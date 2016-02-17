module resids
  ! Residual analysis module

  use log, only: log_log_info
  use numerics, only: numerics_trapz
  use stats, only: stats_mean, stats_variance, stats_fivenum, &
       stats_residuals, stats_mean_sq_err, stats_mean_abs_err, &
       stats_median

  use progvars
  use gaussian, only: gaussian_p
  use output, only: logfile_unit=>output_logfile_unit

  implicit none

  private

  public :: resids_calculate
  public :: resids_stats

contains

  subroutine resids_calculate()

    integer(ip) :: i_p
    real(fp) :: p

    call log_log_info("Constructing residuals.", logfile_unit)

    do i_p = 1, vd_np
       p = vdx%p_range(i_p)
       theor_np_arr(i_p) = abs(gaussian_p(p))**2
    end do

    call stats_residuals(vdx%vd_p_arr, theor_np_arr, resid_np_arr)
    resid_np_mask(:) = resid_np_arr(:) .ge. resid_p_eps

    do i_p = 1, vd_np
       resid_np_cum_arr(i_p) = numerics_trapz(resid_np_arr(1:i_p), vdx%dp)
    end do

  end subroutine resids_calculate

  subroutine resids_stats()
    call log_log_info("Generating statistical summary of residuals.", logfile_unit)

    resid_mean = stats_mean(resid_np_arr, mask=resid_np_mask)
    resid_var = stats_variance(resid_np_arr, mask=resid_np_mask)

    call stats_fivenum(resid_np_arr, resid_fivenum_arr, mask=resid_np_mask)
    resid_mean_sq_err = stats_mean_sq_err(resid_np_arr, mask=resid_np_mask)
    resid_mean_abs_err = stats_mean_abs_err(resid_np_arr, mask=resid_np_mask)

  end subroutine resids_stats

end module resids
