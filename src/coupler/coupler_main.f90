program coupler_main
  use constants_mod,           only: constants_init

  use time_manager_mod,        only: time_type, set_calendar_type, set_time
  use time_manager_mod,        only: set_date, get_date, days_in_month, month_name
  use time_manager_mod,        only: operator(+), operator(-), operator (<)
  use time_manager_mod,        only: operator (>), operator ( /= ), operator ( / )
  use time_manager_mod,        only: operator (*), THIRTY_DAY_MONTHS, JULIAN
  use time_manager_mod,        only: NOLEAP, NO_CALENDAR, INVALID_CALENDAR
  use time_manager_mod,        only: date_to_string, increment_date
  use time_manager_mod,        only: operator(>=), operator(<=), operator(==)

  use fms_mod,                 only: open_namelist_file, field_exist, file_exist, check_nml_error
  use fms_mod,                 only: uppercase, error_mesg, write_version_number
  use fms_mod,                 only: fms_init, fms_end, stdout

  use fms_io_mod,              only: fms_io_exit

  use diag_manager_mod,        only: diag_manager_init, diag_manager_end, diag_grid_end
  use diag_manager_mod,        only: get_base_date
  use diag_manager_mod,        only: diag_manager_set_time_end

  use field_manager_mod,       only: MODEL_ATMOS

  use tracer_manager_mod,      only: tracer_manager_init, get_tracer_index
  use tracer_manager_mod,      only: get_number_tracers, get_tracer_names, NO_TRACER

  use mpp_mod,                 only: mpp_init, mpp_pe, mpp_npes, mpp_root_pe
  use mpp_mod,                 only: stderr, stdlog, mpp_error, NOTE, FATAL, WARNING
  use mpp_mod,                 only: mpp_set_current_pelist, mpp_declare_pelist
  use mpp_mod,                 only: input_nml_file

  use mpp_io_mod,              only: mpp_open, mpp_close, mpp_io_clock_on
  use mpp_io_mod,              only: MPP_NATIVE, MPP_RDONLY, MPP_DELETE

  use atmosphere_mod,          only: atmosphere_init, atmosphere_end,    &
                                     atmosphere_dynamics, atmosphere_state_update

  implicit none

!-----------------------------------------------------------------------

  character(len=128) :: version = '$Id$'
  character(len=128) :: tag = '$Name$'

!-----------------------------------------------------------------------
!---- model defined-types ----

!-----------------------------------------------------------------------
! ----- coupled model time -----

  type (time_type) :: Time, Time_init, Time_end, Time_step
  integer :: num_atmos_calls, na

! ----- coupled model initial date -----

  integer :: date_init(6) = (/ 0, 0, 0, 0, 0, 0 /)
  integer :: calendar_type = INVALID_CALENDAR

!-----------------------------------------------------------------------
!------ namelist interface -------

  integer, dimension(6) :: restart_interval = (/ 0, 0, 0, 0, 0, 0/)
  integer, dimension(6) :: current_date     = (/ 0, 0, 0, 0, 0, 0 /)
  character(len=17) :: calendar = '                 '

  integer :: months=0, days=0, hours=0, minutes=0, seconds=0
  integer :: dt_atmos = 0  ! fluxes passed between atmosphere & ice/land

  integer :: atmos_npes=0

  namelist /coupler_nml/ current_date, calendar,                          &
                         months, days, hours, minutes, seconds, dt_atmos, &
                         atmos_npes

  integer :: id_atmos_model_init

  character(len=48), parameter                    :: mod_name = 'coupler_main_mod'
 
  real :: dsec
  integer :: nc
  character(len=128) :: text

!#######################################################################

  call mpp_init()
  
  call fms_init
  call constants_init
  call coupler_init

  Time = Time_init

  do na = 1, num_atmos_calls

     call atmosphere_dynamics (Time, Time_step)
     call atmosphere_state_update (Time)

     Time = Time + Time_step

  enddo ! end of na (fast loop)

  call coupler_end

  call fms_end

!-----------------------------------------------------------------------

contains

!#######################################################################

  subroutine coupler_init

!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------

    character(len=64), parameter    :: sub_name = 'coupler_init'
    integer :: unit,  ierr, io, m, i, outunit, logunit, errunit
    integer :: date(6)
    type (time_type) :: Run_length
    character(len=9) :: month
    integer :: pe, npes

    integer :: n
    character(len=256) :: err_msg
    character(len=64)  :: filename, fieldname

!-----------------------------------------------------------------------

    outunit = stdout()
    errunit = stderr()
    logunit = stdlog()

!----- write version to logfile -------
    call write_version_number(version, tag)

!----- read namelist -------

    unit = open_namelist_file()
    ierr=1; do while (ierr /= 0)
       read  (unit, nml=coupler_nml, iostat=io, end=10)
       ierr = check_nml_error (io, 'coupler_nml')
    enddo
10  call mpp_close(unit)


        if ( sum(current_date) <= 0 ) then
            call error_mesg ('program coupler',  &
                 'no namelist value for base_date or current_date', FATAL)
        else
            date  = current_date
        endif

!----- override calendar type with namelist value -----

        select case( uppercase(trim(calendar)) )
        case( 'JULIAN' )
            calendar_type = JULIAN
        case( 'NOLEAP' )
            calendar_type = NOLEAP
        case( 'THIRTY_DAY' )
            calendar_type = THIRTY_DAY_MONTHS
        case( 'NO_CALENDAR' )
            calendar_type = NO_CALENDAR
        end select

    call set_calendar_type (calendar_type, err_msg)
    if(err_msg /= '') then
      call mpp_error(FATAL, 'ERROR in coupler_init: '//trim(err_msg))
    endif

!------ initialize diagnostics manager ------

    call diag_manager_init

!------ reset pelist to "full group" ------

    call mpp_set_current_pelist()

!----- set initial and current time types ------

    Time_init = set_date (date(1), date(2), date(3),  &
                          date(4), date(5), date(6))

!----- compute the ending time -----

    Time_end = Time_init
    do m=1,months
       Time_end = Time_end + set_time(0,days_in_month(Time_end))
    end do
    Time_end = Time_end + set_time(hours*3600+minutes*60+seconds, days)

    call diag_manager_set_time_end(Time_end)

    Run_length = Time_end - Time_init

!-----------------------------------------------------------------------
!----- compute the time steps ------

    Time_step = set_time (dt_atmos,0)

!----- determine maximum number of iterations per loop ------

    num_atmos_calls = Run_length / Time_step

!-----------------------------------------------------------------------

    call atmosphere_init (Time_init, Time_init, Time_step, Time_end)

!-----------------------------------------------------------------------

  end subroutine coupler_init

!#######################################################################

  subroutine coupler_end

    call atmosphere_end (Time)

    call fms_io_exit
    call diag_manager_end (Time)
    call mpp_set_current_pelist()

!-----------------------------------------------------------------------

  end subroutine coupler_end

!#######################################################################

  end program coupler_main

