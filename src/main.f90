program main
  ! Main executable

  ! Imports -- library dependencies
  use log, only: log_log_info
  use string, only: string_val
  use stats, only: stats_residuals
  use numerics, only: numerics_trapz

  ! Imports -- program modules
  use progvars
  use setup, only: setup_init, setup_cleanup
  use output, only: output_psi_xt, output_vd_counts, output_vd_residuals, &
       output_vd_t, logfile_unit=>output_logfile_unit
  use propagate, only: propagate_psi
  use gaussian, only: gaussian_p

  implicit none

  integer(ip) :: i_t

  integer(ip) :: i_p
  real(fp) :: p

  call setup_init()

  call log_log_info("Beginning time propagation.", logfile_unit)
  do i_t = 1, nt

     call propagate_psi(psi_arr, i_t)
     call vdx%update(psi_arr)

     if (mod(i_t, print_mod_t) .eq. 0) then
        call log_log_info("Timestep "//string_val(i_t)//" of "// &
             string_val(nt), logfile_unit)
        call output_vd_t()
        call output_psi_xt()
     end if

  end do
  call log_log_info("Time propagation complete.", logfile_unit)

  call vdx%normalize()
  call output_vd_counts()

  call log_log_info("Constructing residuals.", logfile_unit)
  do i_p = 1, vd_np
     p = vdx%p_range(i_p)
     theor_np_arr(i_p) = abs(gaussian_p(p))**2
  end do
  call stats_residuals(vdx%vd_p_arr, theor_np_arr, resid_np_arr)

  do i_p = 1, vd_np
     resid_np_cum_arr(i_p) = numerics_trapz(resid_np_arr(1:i_p), vdx%dp)
  end do
  call output_vd_residuals()

  call setup_cleanup()
end program main
