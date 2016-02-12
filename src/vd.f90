module vd

  use log, only: log_log_critical, log_stderr
  use numerics, only: numerics_cmplx_phase, numerics_d1, numerics_d2, &
       numerics_trapz
  use dists, only: dists_gaussian

  use progvars
  use output, only: output_logfile_unit

  implicit none

  private

  public :: vd_init
  public :: vd_cleanup
  public :: vd_normalize
  public :: vd_update

  integer(ip), parameter :: logfile_unit = output_logfile_unit

  ! Work arrays
  real(fp), allocatable :: phi_arr(:), p_arr(:)
  real(fp), allocatable :: p_var_arr(:)
  real(fp), allocatable :: mag_arr(:), j_arr(:)
  real(fp), allocatable :: prob_arr(:)

contains

  subroutine vd_init()
    logical :: sane

    character(:), allocatable :: error_msg
    
    ! Check internal / external / virtual detector grid limits

    sane = .true.

    if (nx .le. nxl_external + nxr_external) then
       sane = .false.
    else if (nxl_external .le. vd_nxl) then
       sane = .false.
    else if (nxr_external .le. vd_nxr) then
       sane = .false.
    end if

    if (.not. sane) then
       error_msg = "Error in numerical / virtual detector grid "// &
            "configuration; exiting abnormally."
       call log_log_critical(error_msg, logfile_unit)
       call log_log_critical(error_msg, log_stderr)
       call exit(0)
    end if

    ! Calculate edges of VD grid
    vd_xl_min = nxl_external - vd_nxl
    vd_xl_max = nxl_external
    vd_xr_min = nx - nxr_external
    vd_xr_max = vd_xr_min + vd_nxr

    ! Allocate work arrays
    allocate(phi_arr(nx))
    allocate(p_arr(nx))
    allocate(p_var_arr(nx))
    allocate(mag_arr(nx))
    allocate(j_arr(nx))
    allocate(prob_arr(nx))

  end subroutine vd_init

  subroutine vd_cleanup
    ! Deallocate work arrays

    deallocate(prob_arr)
    deallocate(j_arr)
    deallocate(mag_arr)
    deallocate(p_var_arr)
    deallocate(p_arr)
    deallocate(phi_arr)
  end subroutine vd_cleanup

  subroutine vd_normalize(np_arr)
    ! Normalize virtual detector spectrum
    real(fp), intent(inout) :: np_arr(:)
    real(fp) :: np_norm

    np_norm = numerics_trapz(np_arr, vd_dp)
    np_arr(:) = np_arr(:) / np_norm

  end subroutine vd_normalize

  subroutine vd_update(psi_arr, np_arr)
    ! Update virtual detector
    complex(fp), intent(in) :: psi_arr(:)
    real(fp), intent(inout) :: np_arr(:)

    call vd_fill_arrays(psi_arr)
    call vd_calc_pj()
    call vd_bin(np_arr)
    
  end subroutine vd_update

  subroutine vd_bin(np_arr)
    real(fp), intent(inout) :: np_arr(:)

    integer(ip) :: i_x

    do i_x = vd_xl_min, vd_xl_max
       call accumulate_counts(i_x, np_arr)
    end do
    do i_x = vd_xr_min, vd_xr_max
       call accumulate_counts(i_x, np_arr)
    end do

  contains

    subroutine accumulate_counts(i_x, count_arr)
      integer(ip), intent(in) :: i_x
      real(fp), intent(inout) :: count_arr(:)

      real(fp) :: p, p_mu, p_var
      real(fp) :: scale
      integer(ip) :: i_p
      ! mean
      p_mu = p_arr(i_x)
      ! variance
      p_var = p_var_arr(i_x)

      if (vd_semi_classical) then
         ! Semi-classical case uses delta-function histogramming
         ! We mimic this by using a Gaussian with variation vd_dp / 5, so that
         ! roughly 100% of the binned distribution will be located at a single
         ! grid point.
         p_var = vd_dp / 16
      end if

      ! Make Gaussian of variance p_var around p, and populate momentum
      ! distribution using that
      scale = dt * abs(j_arr(i_x))

      do i_p = 1, vd_np
         p = vd_p_range(i_p)
         count_arr(i_p) = count_arr(i_p) + scale * dists_gaussian(p, p_mu, p_var)
      end do

    end subroutine accumulate_counts

  end subroutine vd_bin

  subroutine vd_calc_pj()

    ! Calculate p_local = d(phi) / dx
    call numerics_d1(phi_arr(vd_xl_min - 1 : vd_xl_max + 1), &
         p_arr(vd_xl_min - 1 : vd_xl_max + 1), dx)
    call numerics_d1(phi_arr(vd_xr_min - 1 : vd_xr_max + 1), &
         p_arr(vd_xr_min - 1 : vd_xr_max + 1), dx)

    ! Calculate second moment of local momentum dist.
    ! From Iafrate, G., et al. Journal de Physique 42.C7.10 (1981).
    call numerics_d2(log(mag_arr(vd_xl_min - 1 : vd_xl_max + 1)), &
         p_var_arr(vd_xl_min - 1 : vd_xl_max + 1), dx)
    call numerics_d2(log(mag_arr(vd_xr_min - 1 : vd_xr_max + 1)), &
         p_var_arr(vd_xr_min - 1 : vd_xr_max + 1), dx)

    p_var_arr(:) = - hbar**2 / 4.0_fp * p_var_arr(:)

    ! Calculate j_local = rho * p_local / m
    j_arr(vd_xl_min : vd_xl_max) = mag_arr(vd_xl_min : vd_xl_max) / m * &
         p_arr(vd_xl_min : vd_xl_max)
    j_arr(vd_xr_min : vd_xr_max) = mag_arr(vd_xr_min : vd_xr_max) / m * &
         p_arr(vd_xr_min : vd_xr_max)

  end subroutine vd_calc_pj

  subroutine vd_fill_arrays(psi_arr)
    ! Calculate R^2(x,t), S(x,t) where Psi = R exp(i S / hbar)
    complex(fp), intent(in) :: psi_arr(:)

    complex(fp) :: z
    integer(ip) :: i_x

    ! Left end of VD grid
    do i_x = vd_xl_min - 1, vd_xl_max + 1
       call fill_by_index(i_x)
    end do

    ! Right end of VD grid
    do i_x = vd_xr_min - 1, vd_xr_max + 1
       call fill_by_index(i_x)
    end do

  contains

    subroutine fill_by_index(i_x)
      integer(ip), intent(in) :: i_x

      z = psi_arr(i_x)

      ! multiply by hbar here since Psi ~ exp(i S / hbar)
      phi_arr(i_x) = numerics_cmplx_phase(z) * hbar
      mag_arr(i_x) = abs(z)**2

    end subroutine fill_by_index

  end subroutine vd_fill_arrays

end module vd
