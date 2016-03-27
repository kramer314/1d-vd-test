module vd

  use ieee_arithmetic

  use numerics, only: numerics_cmplx_phase, numerics_d1, numerics_d2
  use log, only: log_log_critical, log_stderr

  use precision, only: ip, fp

  implicit none

  private

  public :: vd_get_local_quantities
  public :: vd_get_indices

contains

  subroutine vd_get_indices(nx, nxl_ext, nxr_ext, nxl_vd, nxr_vd, xl_min, &
       xl_max, xr_min, xr_max)
    ! Get VD indices relative to total spatial grid
    !
    ! This method also checks whether the external and VD region sizes are
    ! coherent, and stops with an error if not.
    !
    ! nx :: spatial grid size
    ! nxl_ext :: number of left external points
    ! nxr_ext :: number of right external points
    ! nxl_vd :: number of left VD points
    ! nxr_vd :: number of right VD points
    ! xl_min :: left VD minimum index (returned value)
    ! xl_max :: left VD maximum index (returned value)
    ! xr_min :: right VD minimum index (returned value)
    ! xr_max :: right VD maximum index (returned value)
    integer, intent(in) :: nx
    integer, intent(in) :: nxl_vd, nxl_ext
    integer, intent(in) :: nxr_vd, nxr_ext
    integer, intent(out) :: xl_min, xl_max
    integer, intent(out) :: xr_min, xr_max

    ! Make sure external and VD region sizes make sense
    call vd_validate_spatial_input(nx, nxl_ext, nxr_ext, nxl_vd, nxr_vd)

    ! We have to be careful here, since Fortran arrays are 1-indexed
    xl_min = nxl_ext + 1
    xl_max = nxl_ext + nxl_vd
    xr_max = nx - nxr_ext
    xr_min = xr_max - (nxr_vd - 1)

  end subroutine vd_get_indices

  subroutine vd_get_local_quantities(psi_arr, dx, m, hbar, p_mu, p_var, j)
    ! Get local p_mu, p_var, j quantities
    !
    ! psi_arr :: 3-element slice of wavefunction surround the virtual detector,
    !   in the desired spatial component direction. We require three elements
    !   because we use three-point finite differencing for derivatives.
    ! dx :: spatial component grid step
    ! m :: mass
    ! hbar :: units
    ! p_mu :: local component momentum value (to be returned)
    ! p_var :: local component momentum variance (to be returned)
    ! j :: local component flux (to be returned)
    complex(fp), intent(in) :: psi_arr(3)
    real(fp), intent(in) :: dx
    real(fp), intent(in) :: m
    real(fp), intent(in) :: hbar
    real(fp), intent(inout) :: p_mu
    real(fp), intent(inout) :: p_var
    real(fp), intent(inout) :: j

    real(fp) :: eps_fp

    integer(ip) :: i_x

    real(fp) :: phi_arr(3), mag_arr(3), p_arr(3), p_var_arr(3)
    complex(fp) :: z

    eps_fp = epsilon(1.0_fp)

    do i_x = 1, 3
       z = psi_arr(i_x)
       phi_arr(i_x) = numerics_cmplx_phase(z) * hbar
       mag_arr(i_x) = abs(z)**2
    end do

    ! Get local (mean) momentum
    call numerics_d1(phi_arr, p_arr, dx)

    p_mu = p_arr(2)

    ! Get local current
    j = mag_arr(2) / m * p_mu

    ! Calculate local momentum variance
    ! Note that this is numerically unstable!
    call numerics_d2(log(mag_arr), p_var_arr, dx)
    p_var = - hbar**2 / 4.0_fp * p_var_arr(2)

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
    if (ieee_is_nan(p_var) .or. (p_var .lt. eps_fp)) then
       p_var = eps_fp
    end if

  end subroutine vd_get_local_quantities

  subroutine vd_validate_spatial_input(nx, nxl_ext, nxr_ext, nxl_vd, nxr_vd)
    ! Validate spatial grid setup and exit abnormally if errors are found
    !
    ! nx :: number of spatial grid points
    ! nxl_ext :: number of left external grid points
    ! nxr_ext :: number of right external grid points
    ! nxl_vd :: number of left virtual detector points
    ! nxr_vd :: number of right virtuald detector points
    integer(ip) :: nx
    integer(ip) :: nxl_ext
    integer(ip) :: nxr_ext
    integer(ip) :: nxl_vd
    integer(ip) :: nxr_vd

    logical :: sane
    character(:), allocatable :: error_msg

    sane = (nxl_ext .gt. 0) .and. (nxr_ext .gt. 0)
    if (.not. sane) then
       error_msg = "No external region defined in numerical grid; "// &
            "stopping abnormally."
       call log_log_critical(error_msg, log_stderr)
       stop 0
    end if

    sane = (nxl_vd .gt. 0) .and. (nxr_vd .gt. 0)
    if (.not. sane) then
       error_msg = "Virtual detectors must be present in numerical grid; "// &
            "stopping abnormally."
       call log_log_critical(error_msg, log_stderr)
       stop 0
    end if

    sane = (nx .gt. nxl_ext + nxr_ext + nxl_vd + nxr_vd)
    if (.not. sane) then
       error_msg = "Numerical grid size greater than number of allocated "// &
            "points; no internal region possible; exiting abnormally."
       call log_log_critical(error_msg, log_stderr)
       stop 0
    end if

  end subroutine vd_validate_spatial_input

end  module vd
