module vd
  ! Virtual detector module

  ! Imports -- intrinsic
  use ieee_arithmetic

  ! Imports -- library dependencies
  use log, only: log_log_critical, log_stderr
  use numerics, only: numerics_cmplx_phase, numerics_d1, numerics_d2, &
       numerics_trapz, numerics_linspace
  use dists, only: dists_gaussian

  ! Imports -- program variables
  use precision, only: ip, fp

  implicit none

  private

  public :: vd_obj
  type vd_obj
     ! Virtual detector object
     logical :: semi_classical

     real(fp) :: dx
     integer(ip) :: nxl, nxr
     integer(ip) :: xl_min, xl_max
     integer(ip) :: xr_min, xr_max

     integer(ip) :: np
     real(fp) :: p_min, p_max
     real(fp) :: dp
     real(fp), allocatable :: p_range(:)
     real(fp), allocatable :: vd_p_arr(:)

     real(fp) :: dt

     real(fp) :: hbar
     real(fp) :: m

     real(fp), allocatable :: phi_arr_l(:), phi_arr_r(:)
     real(fp), allocatable :: p_arr_l(:), p_arr_r(:)
     real(fp), allocatable :: p_var_arr_l(:), p_var_arr_r(:)
     real(fp), allocatable :: mag_arr_l(:), mag_arr_r(:)
     real(fp), allocatable :: j_arr_l(:), j_arr_r(:)

   contains
     procedure :: init => vd_init
     procedure :: cleanup => vd_cleanup
     procedure :: normalize => vd_normalize
     procedure :: update => vd_update
     procedure :: bin => vd_bin
     procedure :: fill_arrays => vd_fill_arrays
  end type vd_obj

contains

  subroutine vd_init(this, nx, nxl_ext, nxr_ext, nxl_vd, nxr_vd, dx, np, &
       p_min, p_max, dt, sc, hbar, m)
    ! Initialize VD object
    class(vd_obj), intent(inout) :: this
    integer(ip), intent(in) :: nx
    integer(ip), intent(in) :: nxl_ext
    integer(ip), intent(in) :: nxr_ext
    integer(ip), intent(in) :: nxl_vd
    integer(ip), intent(in) :: nxr_vd
    real(fp), intent(in) :: dx
    integer(ip), intent(in) :: np
    real(fp), intent(in) :: p_min
    real(fp), intent(in) :: p_max
    real(fp), intent(in) :: dt
    logical, intent(in) :: sc
    real(fp), intent(in) :: hbar
    real(fp), intent(in) :: m

    logical :: sane
    character(:), allocatable :: error_msg

    ! Validate spatial grid parameters
    sane = (nx .gt. nxl_ext + nxr_ext + nxl_vd + nxr_vd)
    if (.not. sane) then
       error_msg = "Error in numerical grid; stopping abnormally"
       call log_log_critical(error_msg, log_stderr)
       deallocate(error_msg)
       stop 0
    end if

    ! Validate momentum grid parameters
    sane = (p_min .lt. p_max)
    if (.not. sane) then
       error_msg = "Error in momentum grid; stopping abnormally"
       call log_log_critical(error_msg, log_stderr)
       deallocate(error_msg)
       stop 0
    end if

    ! Assign VD scalar parameters
    this%dx = dx
    this%nxl = nxl_vd
    this%nxr = nxr_vd

    this%np = np
    this%p_min = p_min
    this%p_min = p_max

    this%dt = dt

    this%semi_classical = sc

    this%hbar = hbar
    this%m = m

    ! Allocate and construct VD counts and momentum grid
    allocate(this%vd_p_arr(np))
    this%vd_p_arr(:) = 0.0_fp
    allocate(this%p_range(np))
    call numerics_linspace(p_min, p_max, this%p_range, this%dp)

    ! Calculate VD grid edge indices on the total spatial grid
    ! Because Fortran arrays are by default 1-indexed, we have to be careful
    ! here.
    this%xl_min = nxl_ext + 1
    this%xl_max = nxl_ext + nxl_vd
    this%xr_max = nx - nxr_ext
    this%xr_min = this%xr_max - nxr_vd + 1

    ! Allocate internal work arrays, including the points on the sides of the
    ! detectors (because we need to compute derivatives).
    ! These are 0-based so that the virtual detector points start at index 1,
    ! which aligns with Fortran's default array indexing.
    allocate(this%phi_arr_l(0 : nxl_vd + 1), this%phi_arr_r(0 : nxr_vd + 1))
    allocate(this%p_arr_l(0 : nxl_vd + 1), this%p_arr_r(0 : nxr_vd + 1))
    allocate(this%p_var_arr_l(0 : nxl_vd + 1), &
         this%p_var_arr_r(0 : nxr_vd + 1))
    allocate(this%mag_arr_l(0 : nxl_vd + 1), this%mag_arr_r(0 : nxr_vd + 1))
    allocate(this%j_arr_l(0 : nxl_vd + 1), this%j_arr_r(0 : nxr_vd + 1))

  end subroutine vd_init

  subroutine vd_cleanup(this)
    ! Deallocate all object arrays
    class(vd_obj), intent(inout) :: this

    deallocate(this%phi_arr_l, this%phi_arr_r)
    deallocate(this%p_arr_l, this%p_arr_r)
    deallocate(this%p_var_arr_l, this%p_var_arr_r)
    deallocate(this%mag_arr_l, this%mag_arr_r)
    deallocate(this%j_arr_l, this%j_arr_r)

    deallocate(this%vd_p_arr)
  end subroutine vd_cleanup

  subroutine vd_normalize(this)
    ! Normalize virtual detector spectrum
    class(vd_obj), intent(inout) :: this

    real(fp) :: np_norm

    np_norm = numerics_trapz(this%vd_p_arr, this%dp)
    this%vd_p_arr(:) = this%vd_p_arr(:) / np_norm
  end subroutine vd_normalize

  subroutine vd_bin(this)
    ! Update virtual detector counts
    class(vd_obj), intent(inout) :: this

    integer(ip) :: i_x_vd

    ! Accumulate counts along left virtual detector grid
    do i_x_vd = 1, this%nxl
       call accumulate_counts(i_x_vd, this%p_arr_l, this%p_var_arr_l, &
            this%j_arr_l)
    end do

    ! Accumulate counts along right virtual detector grid
    do i_x_vd = 1, this%nxr
       call accumulate_counts(i_x_vd, this%p_arr_r, this%p_var_arr_r, &
            this%j_arr_r)
    end do

  contains

    subroutine accumulate_counts(i_x_vd, p_arr, p_var_arr, j_arr)
      ! Method that actually updates virtual detector counts
      !
      ! i_x :: spatial grid index
      ! p_arr :: local momentum array, indexed by i_x
      ! p_var_arr :: local variance array, indexed by i_x
      ! j_arr :: local current array, indexed by i_x
      integer(ip), intent(in) :: i_x_vd
      real(fp), intent(in) :: p_arr(:)
      real(fp), intent(in) :: p_var_arr(:)
      real(fp), intent(in) :: j_arr(:)

      integer(ip) :: i_p
      integer(ip) :: i_p_sc

      real(fp) :: p
      real(fp) :: p_mu
      real(fp) :: p_var

      real(fp) :: scale
      real(fp) :: g

      ! mean momentum
      p_mu = p_arr(i_x_vd)
      ! momentum variance
      p_var = p_var_arr(i_x_vd)

      ! semi-classical adjustments
      ! The semi-classical case uses delta-function histogramming, which we
      ! mimic by using a Gaussian with variance dp/4 so that roughly 100% of
      ! the binned distribution will be located at a single grid point.
      if (this%semi_classical) then
         ! shift p_mu to be at the nearest grid point
         i_p_sc = nint( (p_mu - this%p_min) / this%dp )
         p_mu = this%p_min + this%dp * i_p_sc

         p_var = this%dp / 16
      end if

      scale = this%dt * abs(j_arr(i_x_vd))

      ! Make gaussian of variance p_var around p_mu, and population momentum
      ! distribution using that.
      do i_p = 1, this%np
         p = this%p_range(i_p)
         g = dists_gaussian(p, p_mu, p_var)
         this%vd_p_arr(i_p) = this%vd_p_arr(i_p) + scale * g
      end do

    end subroutine accumulate_counts

  end subroutine vd_bin

  subroutine vd_update(this, psi_arr)
    ! Update virtual detector counts
    !
    ! this :: vd_obj instance
    ! psi_arr :: 1D slice of wavefunction along the same grid as vd_obj
    !
    ! For multi-dimensional propagations, this routine should be called for
    ! each slice along the grid of vd_obj
    class(vd_obj), intent(inout) :: this
    complex(fp), intent(in) :: psi_arr(:)

    call this%fill_arrays(psi_arr)
    call this%bin()
  end subroutine vd_update

  subroutine vd_fill_arrays(this, psi_arr)
    class(vd_obj), intent(inout) :: this
    complex(fp), intent(in) :: psi_arr(:)

    integer(ip) :: i_x
    integer(ip) :: i_x_vd
    complex(fp) :: z
    real(fp) :: eps_fp
    real(fp) :: p_var

    eps_fp = epsilon(1.0_fp)

    ! left side of VD grid
    do i_x = this%xl_min - 1, this%xl_max + 1
       ! i_x indexes the total spatial grid; we need to convert this to the
       ! internal grids we use.
       i_x_vd = i_x - (this%xl_min - 1)

       z = psi_arr(i_x)
       this%phi_arr_l(i_x_vd) = numerics_cmplx_phase(z) * this%hbar
       this%mag_arr_l(i_x_vd) = abs(z)**2
    end do

    call numerics_d1(this%phi_arr_l(:), this%p_arr_l(:), this%dx)
    this%j_arr_l(:) = this%mag_arr_l(:) / this%m * this%p_arr_l(:)

    ! Calculate variance of local momentum distribution.
    ! Note that this is numerically unstable! Right now, we check for this and
    ! manually attempt to correct for it, but we should think about finding
    ! more stable methods to calculate this or find some alternative measure of
    ! variance.
    call numerics_d2(log(this%mag_arr_l(:)), this%p_var_arr_l(:), this%dx)
    this%p_var_arr_l(:) = -this%hbar**2 / 4.0_fp * this%p_var_arr_l(:)

    ! There are cases in which the values of mag_arr are less than machine
    ! precision so that the derivative operation above is completely
    ! numerically unstable. In these cases, we attempt to solve this by setting
    ! the variance equal to the machine epsilon. This can be justified
    ! (informally) using L'Hopital's rule:
    !
    ! d^2/dx^2 log(rho) = 1 / rho^2 ( rho * rho'' - rho''^2)
    !
    ! If rho <= eps locally then rho' and rho'' <= rho, so
    ! d^2/dx^2 log(rho) <= eps, which means we can't have a better estimate
    ! than setting p_var to be the machine epsilon.
    do i_x_vd = 0, this%nxl + 1
       p_var = this%p_var_arr_l(i_x_vd)
       if (ieee_is_nan(p_var) .or. (p_var .lt. eps_fp)) then
          this%p_var_arr_l(i_x_vd) = eps_fp
       end if
    end do

    ! right side of VD grid -- same procedure as the left side above
    do i_x = this%xr_min - 1, this%xr_max + 1
       i_x_vd = i_x - (this%xr_min - 1)

       z = psi_arr(i_x)
       this%phi_arr_r(i_x_vd) = numerics_cmplx_phase(z) * this%hbar
       this%mag_arr_r(i_x_vd) = abs(z)**2
    end do

    call numerics_d1(this%phi_arr_r(:), this%p_arr_r(:), this%dx)
    this%j_arr_r(:) = this%mag_arr_r / this%m * this%p_arr_r(:)

    call numerics_d2(log(this%mag_arr_r(:)), this%p_var_arr_r(:), &
         this%dx)
    this%p_var_arr_r(:) = -this%hbar**2 / 4.0_fp * this%p_var_arr_r(:)

    do i_x_vd = 0, this%nxr + 1
       p_var = this%p_var_arr_r(i_x_vd)
       if (ieee_is_nan(p_var) .or. (p_var .lt. eps_fp)) then
          this%p_var_arr_r(i_x_vd) = eps_fp
       end if
    end do

  end subroutine vd_fill_arrays

end module vd
