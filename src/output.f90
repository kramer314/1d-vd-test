module output

  use files, only: files_ensure_dir
  use log, only: log_log_info

  use progvars

  implicit none

  private

  public :: output_init
  public :: output_cleanup

  public :: output_psi_xt
  public :: output_vd_counts
  
  public :: output_logfile_unit

  ! Output file unit numers
  integer(ip), parameter :: logfile_unit = 99
  integer(ip), parameter :: output_logfile_unit = logfile_unit
  
  integer(ip), parameter :: psi_xt_unit = 98
  integer(ip), parameter :: vd_px_unit = 97

contains
  subroutine output_init()
    
    call files_ensure_dir(output_dir)

    open(unit=logfile_unit, file=trim(output_dir)//trim(log_fname))
    open(unit=psi_xt_unit, file=trim(output_dir)//trim(psi_xt_fname))
    open(unit=vd_px_unit, file=trim(output_dir)//trim(vd_px_fname))

  end subroutine output_init

  subroutine output_cleanup()
    
    close(unit=logfile_unit)
    close(unit=psi_xt_unit)

  end subroutine output_cleanup
  
  subroutine output_psi_xt()
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
    integer(ip) :: i_px
    real(fp) :: px

    call log_log_info("Writing out VD results", logfile_unit)
    do i_px = 1, size(npx_arr)
       px = vd_px_arr(i_px)
       write(vd_px_unit, *) px, npx_arr(i_px)
    end do
  end subroutine output_vd_counts
end module output
