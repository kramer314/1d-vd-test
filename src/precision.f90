module precision

  use globvars, only: dp, dp_format, dp_format_raw, sp, sp_format, &
       sp_format_raw, ip, ip_format, ip_format_raw

  implicit none

  ! Real precision kind parameter / formatting
  integer(ip), parameter :: fp = dp
  character(*), parameter :: fp_format = dp_format
  character(*), parameter :: fp_format_raw = dp_format_raw
  
end module precision
