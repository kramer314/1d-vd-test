module vd_1d
  use numerics, only: numerics_linspace, numerics_trapz
  use dists, only: dists_gaussian
  use precision, only: ip, fp

  use vd, only: vd_get_local_quantities, vd_get_indices

  implicit none

  private

  public vd_1d_obj
  type vd_1d_obj
     ! 1-dimensional VD point

     ! Array of momentum counts; we have to use pointers here because
     ! allocatable arrays can't be used within derived types
     real(fp), pointer :: vd_p_arr(:)

     ! Total probability flux through virtual detector
     real(fp) :: net_flux

     ! Semi-classical binning flag
     logical :: semi_classical

     ! Spatial grid step
     real(fp) :: dx

     ! Momentum grid parameters
     integer(ip) :: np
     real(fp) :: p_min, p_max
     real(fp) :: dp

     ! Pointer to momentum grid (owned by vd_1d_manager); we have to use
     ! pointers here because allocatable arrays can't be used within derived
     ! types
     real(fp), pointer :: p_range(:)

     ! Temporal grid step
     real(fp) :: dt

     ! Units
     real(fp) :: hbar
     real(fp) :: m

   contains
     ! Initialize object
     procedure :: init => vd_1d_obj_init
     ! Cleanup object, deallocating memory
     procedure :: cleanup => vd_1d_obj_cleanup
     ! Update virtual detector counts
     procedure :: update => vd_1d_obj_update

  end type vd_1d_obj

  public vd_1d_manager
  type vd_1d_manager
     ! 1-dimensional VD manager

     ! Virtual detector result; array of momentum counts
     real(fp), allocatable :: vd_p_arr(:)

     ! Semi-classical binning flag
     logical :: semi_classical

     ! Number of left / right virtual detector points
     integer(ip) :: nxl, nxr

     ! Spatial grid parameters
     real(fp) :: dx
     integer(ip) :: nx

     ! Spatial grid index bounds for left / right virtual detecor points
     integer(ip) :: xl_min, xl_max
     integer(ip) :: xr_min, xr_max

     ! Momentum grid parameters
     integer(ip) :: np
     real(fp) :: p_min, p_max
     real(fp) :: dp

     ! Momentum grid; we use an allocatable pointer here since the VD objects
     ! this class managers need this array as well, but we don't want to
     ! allocate a new array for each VD object.
     real(fp), pointer :: p_range(:)

     ! Temporal grid step
     real(fp) :: dt

     ! Units
     real(fp) :: hbar
     real(fp) :: m

     ! Total probability flux through all virtual detectors
     real(fp) :: net_flux

     ! Array of left/right virtual detector point objects
     type(vd_1d_obj), allocatable :: vdl_arr(:), vdr_arr(:)

   contains
     ! Initialize object
     procedure :: init => vd_1d_manager_init
     ! Cleanup object, deallocating memory
     procedure :: cleanup => vd_1d_manager_cleanup
     ! Update virtual detector point counts
     procedure :: update => vd_1d_manager_update
     ! Combine virtual detector point counts and obtain final result
     procedure :: finalize => vd_1d_manager_finalize

  end type vd_1d_manager

contains
  subroutine vd_1d_manager_init(this, nx, nxl_ext, nxr_ext, nxl_vd, nxr_vd, dx, np, &
       p_min, p_max, dt, sc, hbar, m)
    ! Initialize 1D VD manager object
    !
    ! This method is exposed as vd_1d_manager%init
    !
    ! this :: vd_1d_manager instance
    ! nx :: length of spatial grid
    ! nxl_ext :: number of left external spatial points
    ! nxr_ext :: number of right external spatial points
    ! nxl_vd :: number of left VD points
    ! nxr_vd :: number of right VD points
    ! dx :: spatial grid step
    ! np :: length of momentum grid
    ! p_min :: lower bound of momentum grid
    ! p_max :: upper bound of momentum grid
    ! dt :: temporal grid step
    ! sc :: semi-classical binning flag
    ! hbar :: hbar units
    ! m :: particle mass
    class(vd_1d_manager), intent(inout) :: this
    integer(ip), intent(in) :: nx
    integer(ip), intent(in) :: nxl_ext, nxr_ext
    integer(ip), intent(in) :: nxl_vd, nxr_vd
    real(fp), intent(in) :: dx
    integer(ip), intent(in) :: np
    real(fp), intent(in) :: p_min, p_max
    real(fp), intent(in) :: dt
    logical, intent(in) :: sc
    real(fp), intent(in) :: hbar, m

    integer(ip) :: i_x

    this%nx = nx
    this%dx = dx
    this%nxl = nxl_vd
    this%nxr = nxr_vd

    this%np = np
    this%p_min = p_min
    this%p_max = p_max

    this%dt=  dt

    this%semi_classical = sc

    this %hbar = hbar
    this%m = m

    this%net_flux = 0.0_fp

    ! Verify and calculate VD grid edge indices on the total spatial grid
    call vd_get_indices(this%nx, nxl_ext, nxr_ext, this%nxl, this%nxr, &
         this%xl_min, this%xl_max, this%xr_min, this%xr_max)

    ! Construct momentum grid and count arrays
    allocate(this%vd_p_arr(this%np))
    this%vd_p_arr(:) = 0.0_fp

    allocate(this%p_range(this%np))
    call numerics_linspace(this%p_min ,this%p_max, this%p_range, this%dp)

    ! Construct VD arrays
    allocate(this%vdl_arr(this%nxl))
    allocate(this%vdr_arr(this%nxr))

    ! Initialize individual VD objects
    do i_x = 1, this%nxl
       call this%vdl_arr(i_x)%init(this%dx, this%np, this%p_min, this%p_max, &
            this%dp, this%p_range, this%dt, this%semi_classical, &
            this%hbar, this%m)
    end do

    do i_x = 1, this%nxr
       call this%vdr_arr(i_x)%init(this%dx, this%np, this%p_min, this%p_max, &
            this%dp, this%p_range, this%dt, this%semi_classical, &
            this%hbar, this%m)
    end do
  end subroutine vd_1d_manager_init

  subroutine vd_1d_manager_update(this, psi_arr)
    ! Update VD counts
    !
    ! This method is exposed as vd_1d_manager%update
    !
    ! this :: vd_1d_manager instance
    ! psi_arr :: wavefunction
    class(vd_1d_manager), intent(inout) :: this
    complex(fp), intent(in) :: psi_arr(:)

    integer(ip) :: i_x
    integer(ip) :: psi_i_x

    ! Update left VD objects
    do i_x = 1, this%nxl
       psi_i_x = this%xl_min + (i_x - 1)
       call this%vdl_arr(i_x)%update(psi_arr(psi_i_x - 1 : psi_i_x + 1))
    end do

    ! Update right VD objects
    do i_x = 1, this%nxr
       psi_i_x = this%xr_min + (i_x - 1)
       call this%vdr_arr(i_x)%update(psi_arr(psi_i_x - 1 : psi_i_x + 1))
    end do

  end subroutine vd_1d_manager_update

  subroutine vd_1d_manager_cleanup(this)
    ! Cleanup / deallocate internal arrays
    !
    ! This method is exposed as vd_1d_manager%cleanup
    !
    ! this :: vd_1d_manager instance
    class(vd_1d_manager), intent(inout) :: this

    integer(ip) :: i_x

    ! Cleanup left VDs
    do i_x = 1, this%nxl
       call this%vdl_arr(i_x)%cleanup()
    end do
    deallocate(this%vdl_arr)

    ! Cleanup right VDs
    do i_x = 1, this%nxr
       call this%vdr_arr(i_x)%cleanup()
    end do
    deallocate(this%vdr_arr)

    deallocate(this%p_range)
    deallocate(this%vd_p_arr)

  end subroutine vd_1d_manager_cleanup

  subroutine vd_1d_manager_finalize(this)
    ! Combine VD point counts and calculate normalized total momentum
    ! distribution counts
    !
    ! this :: vd_1d_manager instance
    class(vd_1d_manager), intent(inout) :: this

    integer(ip) :: i_x

    real(fp) :: weight_i_x
    real(fp) :: norm

    ! Get total flux
    do i_x = 1, this%nxl
       this%net_flux = this%net_flux + this%vdl_arr(i_x)%net_flux
    end do

    do i_x = 1, this%nxr
       this%net_flux = this%net_flux + this%vdr_arr(i_x)%net_flux
    end do

    ! Combine results
    do i_x = 1, this%nxl
       weight_i_x = this%vdl_arr(i_x)%net_flux / this%net_flux
       this%vd_p_arr(:) = this%vd_p_arr(:) + &
            weight_i_x * this%vdl_arr(i_x)%vd_p_arr(:)
    end do

    do i_x = 1, this%nxr
       weight_i_x = this%vdr_arr(i_x)%net_flux / this%net_flux
       this%vd_p_arr(:) = this%vd_p_arr(:) + &
            weight_i_x * this%vdr_arr(i_x)%vd_p_arr(:)
    end do

    ! Normalize distribution
    norm = numerics_trapz(this%vd_p_arr, this%dp)
    this%vd_p_arr(:) = this%vd_p_arr(:) / norm

  end subroutine vd_1d_manager_finalize

  subroutine vd_1d_obj_init(this, dx, np, p_min, p_max, dp, p_range, dt, sc, &
       hbar, m)
    ! Initialize virtual detector object.
    !
    ! This method is publicly exposed as vd_obj%init
    !
    ! this :: vd_obj instance
    ! dx :: spatial grid step
    ! np :: number of momentum grid points
    ! p_min :: lower momentum grid boundary
    ! p_max :: upper momentum grid boundary
    ! p_min :: pointer array referenced to momentum grid
    ! dt :: temporal grid step
    ! sc :: semi-classical flag
    ! hbar :: units
    ! m :: particle mass
    class(vd_1d_obj), intent(inout) :: this
    real(fp), intent(in) :: dx
    integer(ip), intent(in) :: np
    real(fp), intent(in) :: p_min
    real(fp), intent(in) :: p_max
    real(fp), intent(in) :: dp
    real(fp), pointer, intent(in) :: p_range(:)
    real(fp), intent(in) :: dt
    logical, intent(in) :: sc
    real(fp), intent(in) :: hbar
    real(fp), intent(in) :: m

    this%dx = dx

    this%np = np
    this%p_min = p_min
    this%p_max = p_max
    this%dp = dp

    this%p_range => p_range

    this%dt = dt

    this%semi_classical = sc

    this%hbar = hbar
    this%m = m

    allocate(this%vd_p_arr(this%np))
    this%vd_p_arr(:) = 0.0_fp

    this%net_flux = 0.0_fp
  end subroutine vd_1d_obj_init

  subroutine vd_1d_obj_cleanup(this)
    ! Cleanup / deallocate internal arrays
    !
    ! This method is exposed as vd_1d_obj%cleanup
    !
    ! this :: vd_1d_obj instance
    class(vd_1d_obj), intent(inout) :: this

    ! Since this%p_range points to an array managed by the VD manager, we can't
    ! deallocate it here. Instead, we just dereference the pointer.
    nullify(this%p_range)

    deallocate(this%vd_p_arr)

  end subroutine vd_1d_obj_cleanup

  subroutine vd_1d_obj_update(this, psi_arr)
    ! Update virtual detector counts
    !
    ! This method is exposed as vd_1d_obj%update
    !
    ! this :: vd_1d_obj instance
    ! psi_arr :: 3-element slice of wavefunction values surrounding vd_obj,
    !   in the direction that vd_obj monitors. We require three elements
    !   because we use three-point finite differencing for derivatives.
    class(vd_1d_obj), intent(inout) :: this
    complex(fp), intent(in) :: psi_arr(3)

    real(fp) :: p_mu
    real(fp) :: p_var
    real(fp) :: j

    integer(ip) :: i_p
    real(fp) :: p
    real(fp) :: scale
    real(fp) :: g

    ! Get local p_mu, p_var, j
    call vd_get_local_quantities(psi_arr, this%dx, this%m, this%hbar, p_mu, &
         p_var, j)

    scale = this%dt * abs(j)

    if (this%semi_classical) then
       ! Shift p_mu to be at the nearest grid point
       i_p = nint( (p_mu - this%p_min) / this%dp )

       if (i_p .ge. 1 .and. i_p .le. this%np) then
          ! Divide scale by dp here so that delta-distribution is normalized
          this%vd_p_arr(i_p) = this%vd_p_arr(i_p) + scale / this%dp
       end if
    else
       ! Make Gaussian of variance p_var around p_mu and populate momentum
       ! distribution using that
       do i_p = 1, this%np
          p = this%p_range(i_p)
          g = dists_gaussian(p, p_mu, p_var)
          this%vd_p_arr(i_p) = this%vd_p_arr(i_p) + scale * g
       end do
    end if

    this%net_flux = this%net_flux + scale

  end subroutine vd_1d_obj_update

end module vd_1d
