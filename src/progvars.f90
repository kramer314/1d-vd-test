module progvars
  use ieee_arithmetic

  use globvars, only: dp, dp_format, dp_format_raw, ip, pi_dp, e_dp, j_dp
  use config, only: config_get_param
  use numerics, only: numerics_linspace

  implicit none

  ! Real kind parameter / formatting
  integer(ip), parameter :: fp = dp
  character(*), parameter :: fp_format = dp_format
  character(*), parameter :: fp_format_raw = dp_format_raw

  ! Numerical constants
  real(fp), parameter :: pi = real(pi_dp, kind=fp)
  real(fp), parameter :: sqrt_pi = sqrt(pi)
  real(fp), parameter :: e = real(e_dp, kind=fp)
  complex(fp), parameter :: j = cmplx(j_dp, kind=fp)

  ! Units
  real(fp) :: hbar

  ! Particle parameters
  real(fp) :: m

  ! Gaussian parameters
  real(fp) :: p0, sig_x, x0

  ! Grid parameters
  real(fp) :: x_min, x_max, dx
  integer(ip) :: nx

  ! Virtual detector parameters
  ! left/right number of grid points outside region of interested
  integer(ip) :: nxl_external, nxr_external
  ! left/right number of virtual detector grid points in external grid
  integer(ip) :: vd_nxl, vd_nxr

  ! semi-classical VD switch
  logical :: vd_semi_classical
  
  ! Useful VD grid indices
  integer(ip) :: vd_xl_min, vd_xl_max, vd_xr_min, vd_xr_max

  ! VD momentum bin parameters
  real(fp) :: vd_p_min, vd_p_max, vd_dp
  integer(ip) :: vd_np

  ! Time grid parameters
  real(fp) :: t_min, t_max, dt
  integer(ip) :: nt

  ! Output parameters
  integer(ip) :: print_mod_x, print_mod_t
  character(:), allocatable :: output_dir
  character(:), allocatable :: log_fname
  character(:), allocatable :: psi_xt_fname
  character(:), allocatable :: vd_p_fname

  ! Residual analysis parameters
  real(fp) :: resid_p_eps

  ! Arrays
  real(fp), allocatable :: x_range(:), t_range(:)
  real(fp), allocatable :: vd_p_range(:)
  real(fp), allocatable :: vd_np_arr(:)
  complex(fp), allocatable :: psi_arr(:)

  real(fp), allocatable :: theor_np_arr(:), resid_np_arr(:)
  real(fp), allocatable :: resid_np_cum_arr(:)

contains

  subroutine progvars_init()
    call progvars_read_params()
    call progvars_allocate_arrays()
    call progvars_set_arrays()
  end subroutine progvars_init

  subroutine progvars_cleanup()
    call progvars_deallocate_arrays()
    call progvars_deallocate_params()
  end subroutine progvars_cleanup

  subroutine progvars_allocate_arrays()
    allocate(psi_arr(nx))
    allocate(t_range(nt))
    allocate(x_range(nx))

    allocate(vd_p_range(vd_np))
    allocate(vd_np_arr(vd_np))

    allocate(theor_np_arr(vd_np))
    allocate(resid_np_arr(vd_np))
    allocate(resid_np_cum_arr(vd_np))
  end subroutine progvars_allocate_arrays

  subroutine progvars_deallocate_arrays()
    deallocate(vd_p_range)
    deallocate(vd_np_arr)

    deallocate(x_range)
    deallocate(t_range)
    deallocate(psi_arr)

    deallocate(theor_np_arr)
    deallocate(resid_np_arr)
    deallocate(resid_np_cum_arr)
  end subroutine progvars_deallocate_arrays

  subroutine progvars_set_arrays()
    ! Initialize numerical grids
    call numerics_linspace(t_min, t_max, t_range, dt)
    call numerics_linspace(x_min, x_max, x_range, dx)

    ! Initialize virtual detector grids
    call numerics_linspace(vd_p_min, vd_p_max, vd_p_range, vd_dp)

    ! Initialize virtual detector counts
    vd_np_arr(:) = 0.0_fp
  end subroutine progvars_set_arrays

  subroutine progvars_read_params()
    ! Read in parameters from input file
    logical :: success

    call config_get_param("hbar", hbar, success)

    call config_get_param("m", m, success)

    call config_get_param("p0", p0, success)
    call config_get_param("x0", x0, success)
    call config_get_param("sig_x", sig_x, success)

    call config_get_param("x_min", x_min, success)
    call config_get_param("x_max", x_max, success)
    call config_get_param("nx", nx, success)

    call config_get_param("t_min", t_min, success)
    call config_get_param("t_max", t_max, success)
    call config_get_param("nt", nt, success)

    call config_get_param("nxl_external", nxl_external, success)
    call config_get_param("nxr_external", nxr_external, success)
    call config_get_param("vd_nxl", vd_nxl, success)
    call config_get_param("vd_nxr", vd_nxr, success)

    call config_get_param("vd_p_min", vd_p_min, success)
    call config_get_param("vd_p_max", vd_p_max, success)
    call config_get_param("vd_np", vd_np, success)

    call config_get_param("vd_semi_classical", vd_semi_classical, success)

    call config_get_param("output_dir", output_dir, success)
    call config_get_param("log_fname", log_fname, success)
    call config_get_param("psi_xt_fname", psi_xt_fname, success)
    call config_get_param("vd_p_fname", vd_p_fname, success)

    call config_get_param("print_mod_x", print_mod_x, success)
    call config_get_param("print_mod_t", print_mod_t, success)

  end subroutine progvars_read_params

  subroutine progvars_deallocate_params()
    deallocate(vd_p_fname)
    deallocate(psi_xt_fname)
    deallocate(log_fname)
    deallocate(output_dir)
  end subroutine progvars_deallocate_params

end module progvars
