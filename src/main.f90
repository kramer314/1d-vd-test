program main
  use log, only: log_log_info
  use string, only: string_val

  use progvars
  use setup, only: setup_init, setup_cleanup
  use output, only: output_psi_xt, output_vd_counts, output_logfile_unit
  use propagate, only: propagate_psi
  use vd, only: vd_update, vd_normalize

  integer(ip) :: logfile_unit = output_logfile_unit
  integer(ip) :: i_t

  call setup_init()

  call log_log_info("Beginning time propagation.", logfile_unit)
  do i_t = 1, nt

     call propagate_psi(psi_arr, i_t)
     call vd_update(psi_arr, npx_arr)

     if (mod(i_t, print_mod_t) .eq. 0) then
        call log_log_info("Timestep "//string_val(i_t)//" of "// &
             string_val(nt), logfile_unit)
        call output_psi_xt()
     end if

  end do
  call log_log_info("Time propagation complete.", logfile_unit)

  call vd_normalize(npx_arr)

  call output_vd_counts()

  call setup_cleanup()
end program main
