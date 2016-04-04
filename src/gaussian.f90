! Copyright (c) 2016 Alex Kramer <kramer.alex.kramer@gmail.com>
! See the LICENSE.txt file at the top-level of this distribution.

module gaussian
  ! 1D Gaussian module

  use progvars

  implicit none

  public :: gaussian_init
  public :: gaussian_cleanup
  public :: gaussian_p
  public :: gaussian_xt

  private

  ! Useful constants
  real(fp) :: p0_m, p0_2m
  real(fp) :: sig_x2
  real(fp) :: sig_p, sig_p2
  real(fp) :: norm_p
  complex(fp) :: j_hb, jhb_m

contains
  subroutine gaussian_init()
    ! Precalculate useful constants
    p0_m = p0 / m
    p0_2m = p0_m / 2.0_fp

    sig_x2 = sig_x ** 2

    j_hb = j / hbar
    jhb_m = j * hbar / m

    sig_p = hbar / sig_x
    sig_p2 = sig_p**2

    norm_p = sqrt(1.0_fp / (sig_p * sqrt_pi))

  end subroutine gaussian_init

  subroutine gaussian_cleanup()
  end subroutine gaussian_cleanup

  real(fp) function gaussian_p(p) result(val)
    ! Taken from Sakurai, p. 58
    real(fp), intent(in) :: p

    real(fp) :: exp_p

    exp_p = exp(- 0.5_fp * ( (p - p0) / sig_p )**2 )

    val = norm_p * exp_p
  end function gaussian_p

  complex(fp) function gaussian_xt(x, t) result(val)
    real(fp), intent(in) :: x
    real(fp), intent(in) :: t

    complex(fp) :: norm_x
    complex(fp) :: exp1, exp2

    norm_x = 1.0_fp / sqrt( sqrt_pi * (sig_x + jhb_m / sig_x * t) )
    exp1 = - ( (x - x0) - p0_m * t )**2 / &
         ( 2.0_fp * sig_x2 * (1 + jhb_m / sig_x2 * t) )
    exp2 = j_hb * p0 * ( (x - x0) - p0_2m * t )

    val = norm_x * exp(exp1 + exp2)

  end function gaussian_xt

end module gaussian
