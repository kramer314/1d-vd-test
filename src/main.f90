program main

  use log, only: log_log_info
  use string, only: string_val
  use stats, only: stats_residuals
  use numerics, only: numerics_trapz

  use progvars
  use setup, only: setup_init, setup_cleanup
  use output, only: output_psi_xt, output_vd_counts, output_vd_t, &
       output_logfile_unit
  use propagate, only: propagate_psi
  use gaussian, only: gaussian_p
  use vd, only: vd_update, vd_normalize

  implicit none

  integer(ip) :: logfile_unit = output_logfile_unit
  integer(ip) :: i_t

  integer(ip) :: i_p
  real(fp) :: p

  call setup_init()

  call log_log_info("Beginning time propagation.", logfile_unit)
  do i_t = 1, nt

     call propagate_psi(psi_arr, i_t)
     call vd_update(psi_arr, vd_np_arr)

     if (mod(i_t, print_mod_t) .eq. 0) then
        call log_log_info("Timestep "//string_val(i_t)//" of "// &
             string_val(nt), logfile_unit)
        call output_vd_t()
        ! call output_psi_xt()
     end if

  end do
  call log_log_info("Time propagation complete.", logfile_unit)

  call vd_normalize(vd_np_arr)

  call log_log_info("Constructing residuals.", logfile_unit)
  do i_p = 1, vd_np
     p = vd_p_range(i_p)
     theor_np_arr(i_p) = abs(gaussian_p(p))**2
  end do
  call stats_residuals(vd_np_arr, theor_np_arr, resid_np_arr)

  do i_p = 1, vd_np
     resid_np_cum_arr(i_p) = numerics_trapz(resid_np_arr(1:i_p), vd_dp)
  end do

  call output_vd_counts()

  call setup_cleanup()
end program main
