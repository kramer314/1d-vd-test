! Copyright (c) 2016 Alex Kramer <kramer.alex.kramer@gmail.com>
! See the LICENSE.txt file at the top-level of this distribution.

module propagate
  use progvars
  use gaussian, only: gaussian_xt, gaussian_init, gaussian_cleanup

  implicit none

  private

  public :: propagate_psi
  public :: propagate_init
  public :: propagate_cleanup

contains

  subroutine propagate_init()
    call gaussian_init()
  end subroutine propagate_init

  subroutine propagate_cleanup()
    call gaussian_cleanup()
  end subroutine propagate_cleanup

  subroutine propagate_psi(psi_arr, i_t)

    complex(fp), intent(inout) :: psi_arr(:)
    integer, intent(in) :: i_t

    real(fp) :: x, t
    integer(ip) :: i_x

    t = t_range(i_t)

    ! For testing efficiency, only fill VD region
    !$omp parallel do private(i_x, x)
    do i_x = vd_xl_min - 1, vd_xl_max + 1
       x = x_range(i_x)
       psi_arr(i_x) = gaussian_xt(x, t)
    end do
    !$omp end parallel do

    !$omp parallel do private(i_x, x)
    do i_x = vd_xr_min - 1, vd_xr_max + 1
       x = x_range(i_x)
       psi_arr(i_x) = gaussian_xt(x, t)
    end do
    !$omp end parallel do

  end subroutine propagate_psi
end module propagate
