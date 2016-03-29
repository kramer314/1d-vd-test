module vd_1d
  use numerics, only: numerics_linspace_index
  use dists, only: dists_gaussian
  use precision, only: ip, fp

  use vd, only: vd_get_local_quantities, vd_validate_quantum_update

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
     ! Number of standard deviations to include in quantum VD binning
     integer(ip) :: vd_np_stdev

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

contains

  subroutine vd_1d_obj_init(this, dx, np, p_min, p_max, dp, p_range, dt, sc, &
       hbar, m, vd_np_stdev)
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
    integer(ip), intent(in), optional :: vd_np_stdev

    this%dx = dx

    this%np = np
    this%p_min = p_min
    this%p_max = p_max
    this%dp = dp

    this%p_range => p_range

    this%dt = dt

    this%semi_classical = sc

    if (present(vd_np_stdev)) then
       this%vd_np_stdev= vd_np_stdev
    end if

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

    real(fp) :: p_min, p_max
    real(fp) :: p_stdev
    integer(ip) :: i_p_min, i_p_max
    logical :: valid

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
       i_p = numerics_linspace_index(p_mu, this%p_range)
       if (i_p .gt. 0_ip) then
          ! Divide scale by dp here so that delta-distribution is normalized
          this%vd_p_arr(i_p) = this%vd_p_arr(i_p) + scale / this%dp
       end if

    else

       ! Populate momentum distribution around p_mu
       p_stdev = sqrt(p_var)
       p_min = p_mu - this%vd_np_stdev * p_stdev
       p_max = p_mu + this%vd_np_stdev * p_stdev

       i_p_min = numerics_linspace_index(p_min, this%p_range)
       i_p_max = numerics_linspace_index(p_max, this%p_range)

       call vd_validate_quantum_update(i_p_min, i_p_max, this%np, valid)

       if (valid) then
          ! Make Gaussian of variance p_var around p_mu and populate momentum
          ! distribution using that
          do i_p = i_p_min, i_p_max
             p = this%p_range(i_p)
             g = dists_gaussian(p, p_mu, p_var)
             this%vd_p_arr(i_p) = this%vd_p_arr(i_p) + scale * g
          end do

       end if

    end if

    this%net_flux = this%net_flux + scale

  end subroutine vd_1d_obj_update

end module vd_1d
