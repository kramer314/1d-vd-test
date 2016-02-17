module setup
  use config, only: config_init, config_cleanup
  use log, only: log_log_info, log_stdout

  use progvars
  use output, only: output_init, output_cleanup, output_logfile_unit
  use propagate, only: propagate_init, propagate_cleanup

  implicit none

  private

  public :: setup_init
  public :: setup_cleanup

  integer(ip), parameter :: logfile_unit = output_logfile_unit

contains
  subroutine setup_init()

    character(120) :: input_fname
    character(:), allocatable :: cmd_to_exec

    ! Get input / config file from command-line argument and read in parameters
    call get_command_argument(1, input_fname)
    call log_log_info("Reading config/input file at: "//trim(input_fname), &
         log_stdout)
    call config_init(trim(input_fname))
    call progvars_init()

    ! Copy input file to output directory
    cmd_to_exec = "cp "//trim(input_fname)//" "//trim(output_dir)
    call log_log_info("Copying config/input file to output directory: "// &
         cmd_to_exec, log_stdout)
    call execute_command_line(cmd_to_exec)

    call log_log_info("Setting up output files. Subsequent logging will "// &
         "be redirected to the specified log file.", log_stdout)
    call output_init()

    call propagate_init()

    call log_log_info("Initialization complete.", logfile_unit)

  end subroutine setup_init

  subroutine setup_cleanup()

    call log_log_info("Beginning module cleanup / exiting.", logfile_unit)

    call propagate_cleanup()
    call output_cleanup()
    call progvars_cleanup()
    call config_cleanup()

    call log_log_info("Cleanup finished; exiting.", log_stdout)
  end subroutine setup_cleanup

end module setup
