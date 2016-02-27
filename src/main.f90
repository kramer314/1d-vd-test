program main
  ! Main executable

  ! Imports -- library dependencies
  use log, only: log_log_info
  use string, only: string_val

  ! Imports -- program modules
  use progvars
  use setup, only: setup_init, setup_cleanup
  use output, only: output_psi_xt, output_vd_counts, output_vd_residuals, &
       output_vd_residual_analysis, output_vd_pt, &
       logfile_unit=>output_logfile_unit
  use propagate, only: propagate_psi
  use resids, only: resids_calculate, resids_stats

  implicit none

  integer(ip) :: i_t

  call setup_init()

  call log_log_info("Beginning time propagation.", logfile_unit)
  do i_t = 1, nt

     call propagate_psi(psi_arr, i_t)
     call vdx%update(psi_arr)

     if (mod(i_t, print_mod_t) .eq. 0) then
        call log_log_info("Timestep "//string_val(i_t)//" of "// &
             string_val(nt), logfile_unit)

        if (write_out_vd_pt) then
           call output_vd_pt()
        end if

        if (write_out_psi_xt) then
           call output_psi_xt()
        end if

     end if

  end do
  call log_log_info("Time propagation complete.", logfile_unit)

  call output_vd_counts()

  call resids_calculate()
  call resids_stats()

  call output_vd_residuals()
  call output_vd_residual_analysis()

  call setup_cleanup()
end program main
