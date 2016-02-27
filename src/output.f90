module output
  ! Output module

  ! Imports -- library dependencies
  use log, only: log_log_info, log_log_critical, log_stderr

  ! Imports -- program modules
  use progvars
  use gaussian, only: gaussian_p

  implicit none

  private

  public :: output_init
  public :: output_cleanup

  public :: output_psi_xt
  public :: output_vd_counts
  public :: output_vd_residuals
  public :: output_vd_residual_analysis
  public :: output_vd_pt

  public :: output_logfile_unit

  ! Output file unit numbers
  integer(ip), parameter :: logfile_unit = 99
  integer(ip), parameter :: output_logfile_unit = logfile_unit

  integer(ip), parameter :: psi_xt_unit = 98
  integer(ip), parameter :: vd_p_unit = 97
  integer(ip), parameter :: vd_resid_unit = 96
  integer(ip), parameter :: vd_resid_analysis_unit = 95
  integer(ip), parameter :: vd_pt_unit = 94


contains
  subroutine output_init()
    ! Initialize output module, including opening output files

    open(unit=logfile_unit, file=trim(output_dir)//trim(log_fname))
    open(unit=vd_p_unit, file=trim(output_dir)//trim(vd_p_fname))
    open(unit=vd_resid_unit, file=trim(output_dir)//trim(vd_resid_fname))
    open(unit=vd_resid_analysis_unit, file=trim(output_dir)//&
         trim(vd_resid_analysis_fname))

    if (write_out_vd_pt) then
       open(unit=vd_pt_unit, file=trim(output_dir)//trim(vd_pt_fname))
    end if

    if (write_out_psi_xt) then
       open(unit=psi_xt_unit, file=trim(output_dir)//trim(psi_xt_fname))
    end if

  end subroutine output_init

  subroutine output_cleanup()
    ! Cleanup output module, including closing output files

    close(unit=logfile_unit)
    close(unit=psi_xt_unit)
    close(unit=vd_p_unit)
    close(unit=vd_pt_unit)
    close(unit=vd_resid_unit)
    close(unit=vd_resid_analysis_unit)

  end subroutine output_cleanup

  subroutine output_psi_xt()
    ! Output time-dependent wavefunction psi(x,t)
    integer(ip) :: i_x
    character(:), allocatable :: error_msg

    if (.not. write_out_psi_xt) then
       error_msg = "write_out_psi_xt=False, but output_psi_xt() called; "// &
            "exiting abnormally."
       call log_log_critical(error_msg, log_stderr)
    end if

    call log_log_info("Writing out psi(x,t)", logfile_unit)
    do i_x = 1, nx
       if (mod(i_x, print_mod_x) .eq. 0) then
          write(psi_xt_unit, fp_format, advance="no") abs(psi_arr(i_x))**2
       end if
    end do
    write(psi_xt_unit, *)

  end subroutine output_psi_xt

  subroutine output_vd_counts()
    ! Output virtual detector momentum-binning counts
    integer(ip) :: i_p
    real(fp) :: p

    call log_log_info("Writing out VD counts", logfile_unit)

    do i_p = 1, size(vdx%vd_p_arr)
       p = vdx%p_range(i_p)
       write(vd_p_unit, "(2"//fp_format_raw//")") p, vdx%vd_p_arr(i_p)
    end do

  end subroutine output_vd_counts

  subroutine output_vd_residuals()
    ! Output virtual detector residual comparisons
    integer(ip) :: i_p
    real(fp) :: p

    call log_log_info("Writing out VD residuals", logfile_unit)
    do i_p = 1, size(vdx%p_range)
       p = vdx%p_range(i_p)
       write(vd_resid_unit, "(4"//fp_format_raw//")") p, theor_np_arr(i_p), &
            resid_np_arr(i_p), resid_np_cum_arr(i_p)
    end do
  end subroutine output_vd_residuals

  subroutine output_vd_residual_analysis()
    ! Output residual statistical analysis

    integer(ip), parameter :: out_unit = vd_resid_analysis_unit

    character(*), parameter :: char_format = "(A)"

    write(out_unit, char_format) "Virtual detector residual statistical analysis"
    write(out_unit, char_format) "Format: [Label] [value]"
    write(out_unit, *)

    write(out_unit, char_format, advance="no") "Cumulative residuals: "
    write(out_unit, fp_format) resid_np_cum_arr(vd_np)

    write(out_unit, char_format, advance="no") "Residual threshold: "
    write(out_unit, fp_format) resid_p_eps

    write(out_unit, char_format, advance="no") "Number of above-threshold residuals "// &
         "used in analysis: "
    write(out_unit, ip_format) count(resid_np_mask)

    write(out_unit, char_format, advance="no") "Min residual value: "
    write(out_unit, fp_format) resid_fivenum_arr(1)

    write(out_unit, char_format, advance="no") "Max residual value: "
    write(out_unit, fp_format) resid_fivenum_arr(5)

    write(out_unit, char_format, advance="no") "Mean residual value: "
    write(out_unit, fp_format) resid_mean

    write(out_unit, char_format, advance="no") "Median residual value: "
    write(out_unit, fp_format) resid_fivenum_arr(3)

    write(out_unit, char_format, advance="no") "Residual variance: "
    write(out_unit, fp_format) resid_var

    write(out_unit, char_format, advance="no") "Lower quartile residual value: "
    write(out_unit, fp_format) resid_fivenum_arr(2)

    write(out_unit, char_format, advance="no") "Upper quartile residual value: "
    write(out_unit, fp_format) resid_fivenum_arr(4)

    write(out_unit, char_format, advance="no") "Mean absolute residual error: "
    write(out_unit, fp_format) resid_mean_abs_err

    write(out_unit, char_format, advance="no") "Mean squared residual error: "
    write(out_unit, fp_format) resid_mean_sq_err

    write(out_unit, *)

  end subroutine output_vd_residual_analysis

  subroutine output_vd_pt()
    ! Output time-dependent virual detector counts
    integer(ip) :: i_p
    real(fp) :: p
    character(:), allocatable :: error_msg

    if (.not. write_out_vd_pt) then
       error_msg = "write_out_vd_pt=False, but output_vd_pt() called; "// &
            "exiting abnormally."
       call log_log_critical(error_msg, log_stderr)
    end if

    call log_log_info("Writing out phi(p; t)", logfile_unit)
    do i_p = 1, size(vdx%p_range)
       p = vdx%p_range(i_p)
       write(vd_pt_unit, fp_format, advance="no") vdx%vd_p_arr(i_p)
    end do
    write(vd_pt_unit, *)

  end subroutine output_vd_pt

end module output
