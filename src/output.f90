module output
  ! Output module

  ! Imports -- library dependencies
  use files, only: files_ensure_dir
  use log, only: log_log_info

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
  public :: output_vd_t

  public :: output_logfile_unit

  ! Output file unit numers
  integer(ip), parameter :: logfile_unit = 99
  integer(ip), parameter :: output_logfile_unit = logfile_unit

  integer(ip), parameter :: psi_xt_unit = 98
  integer(ip), parameter :: vd_p_unit = 97
  integer(ip), parameter :: vd_resid_unit = 96
  integer(ip), parameter :: vd_pt_unit = 95

contains
  subroutine output_init()
    ! Initialize output module, including opening output files

    call files_ensure_dir(output_dir)

    open(unit=logfile_unit, file=trim(output_dir)//trim(log_fname))
    open(unit=psi_xt_unit, file=trim(output_dir)//trim(psi_xt_fname))
    open(unit=vd_p_unit, file=trim(output_dir)//trim(vd_p_fname))
    open(unit=vd_resid_unit, file=trim(output_dir)//trim(vd_resid_fname))
    open(unit=vd_pt_unit, file=trim(output_dir)//trim(vd_pt_fname))

  end subroutine output_init

  subroutine output_cleanup()
    ! Cleanup output module, including closing output files

    close(unit=logfile_unit)
    close(unit=psi_xt_unit)
    close(unit=vd_p_unit)
    close(unit=vd_pt_unit)
    close(unit=vd_resid_unit)

  end subroutine output_cleanup

  subroutine output_psi_xt()
    ! Output time-dependent wavefunction psi(x,t)
    integer(ip) :: i_x

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

  subroutine output_vd_t()
    ! Output time-dependent virual detector counts
    integer(ip) :: i_p
    real(fp) :: p

    call log_log_info("Writing out phi(p; t)", logfile_unit)
    do i_p = 1, size(vdx%p_range)
       p = vdx%p_range(i_p)
       write(vd_pt_unit, fp_format, advance="no") vdx%vd_p_arr(i_p)
    end do
    write(vd_pt_unit, *)

  end subroutine output_vd_t

end module output
