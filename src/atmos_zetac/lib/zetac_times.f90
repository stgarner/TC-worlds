module zetac_times_mod

!-----------------------------------------------------------------------
! module information
!-----------------------------------------------------------------------

use fms_mod,          only : file_exist, check_nml_error,              &
                             mpp_pe, mpp_root_pe, stdout,              &
                             stdlog, write_version_number,             &
                             open_namelist_file,                       &
                             open_file, close_file,                    &
                             error_mesg, FATAL
use mpp_mod,          only : input_nml_file

use time_manager_mod, only : time_type, set_calendar_type,             &
                             operator(+), operator(-),                 &
                             operator(<), operator(>),                 &
                             get_time, set_time, get_date, set_date,   &
                             days_in_month, month_name,                &
                             thirty_day_months, julian,                &
                             NOLEAP, no_calendar

!-----------------------------------------------------------------------
! public interfaces
!-----------------------------------------------------------------------

implicit none
private
public get_times, write_timestamps

logical :: do_netcdf_restart, force_date_from_namelist
integer :: base_year=0
integer :: current_date(6)=0

!-----------------------------------------------------------------------
! version control information
!-----------------------------------------------------------------------

character(len=128)  :: version =  '$Id: zetac_times.f90,v 1.1.2.4.2.5 2005/01/29 22:35:44 stg Exp $'
character(len=128)  :: tag     =  '$Name:  $'

contains

!#######################################################################

subroutine get_times (                                                 &
                                                  Time_init, Time_end, &
                                                           seg, nsegs )
    
type (time_type),                                  intent (out) ::     &
                                                  Time_init, Time_end

integer, optional,                                 intent (in)  ::     &
                                                           seg, nsegs

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

character(len=9)  :: mname
character(len=17) :: calendar

integer :: months=0, days=0, hours=0, minutes=0, seconds=0
integer :: dt_atmos
logical :: read_data=.false., use_data=.false.

namelist /coupler_nml/                                                 &
  current_date,                                                        &
  months, days, hours, minutes, seconds,                               &
  dt_atmos,                                                            &
  calendar

integer :: run_time_days_beg, run_time_secs_beg
integer :: run_time_days_end, run_time_secs_end
real    :: run_time, run_time_beg, run_time_end

integer :: calendar_type=-99
integer :: ierr, unit, io
integer :: num_cpld, num_ocean, num_atmos
integer :: m, nsecs, ndays

integer :: segment, num_segs

!-----------------------------------------------------------------------
! define information for segments
!-----------------------------------------------------------------------
 
  if ( present(seg) ) then
    segment = seg
    num_segs = nsegs
  else
    segment = 1
    num_segs = 1
  endif
 
!-----------------------------------------------------------------------
! read namelist file
!-----------------------------------------------------------------------

   read (input_nml_file, nml=coupler_nml, iostat=io)
   ierr = check_nml_error(io,'coupler_nml')

!-----------------------------------------------------------------------
! write version number and namelists to logfile
!-----------------------------------------------------------------------

  call write_version_number (version, tag)
  if (mpp_pe() == mpp_root_pe()) write (stdlog(), nml=coupler_nml )

!-----------------------------------------------------------------------
! define calendar (unless base_year = 0)
!-----------------------------------------------------------------------

  base_year = current_date(1)
  
  if ( base_year > 0 ) then

     if (calendar(1:6) == 'julian') then
        calendar_type = julian
     else if (calendar(1:6) == 'NOLEAP') then
        calendar_type = NOLEAP
     else if (calendar(1:10) == 'thirty_day') then
        calendar_type = thirty_day_months
     else if (calendar(1:11) == 'no_calendar') then
        calendar_type = no_calendar
     else 
        call error_handler ( 'invalid namelist value for calendar' )
     endif

     call set_calendar_type (calendar_type)
    
!-----------------------------------------------------------------------
! set initial time
!-----------------------------------------------------------------------

     Time_init = set_date ( current_date(1), current_date(2),          &
                            current_date(3), current_date(4),          &
                            current_date(5), current_date(6) )
     if ( current_date(2) < 1 .or. current_date(3) < 1 ) then
        call error_handler ( "initial month and/or day is invalid" )
     endif

  else

     nsecs = (current_date(4)*60 + current_date(5))*60 + current_date(6)
     ndays =  current_date(3)
     Time_init  = set_time (nsecs, ndays)

  endif

  if ( current_date(4) >= 24 ) then
     call error_handler ( "initial hour is invalid" )
  endif

  call get_time ( Time_init, nsecs, ndays )
  run_time_beg = real(nsecs) + 86400.*ndays

!-----------------------------------------------------------------------
! set final time and check for errors
!-----------------------------------------------------------------------

  Time_end = Time_init
  do m=1,months
     Time_end = Time_end + set_time ( 0, days_in_month(Time_end) )
  enddo
  Time_end = Time_end + set_time ( (hours*60+minutes)*60+seconds, days )

  call get_time ( Time_end, nsecs, ndays )
  run_time_end = real(nsecs) + 86400.*ndays

!-----------------------------------------------------------------------
! compute beginning and ending times of segment
!-----------------------------------------------------------------------

  run_time = (run_time_end - run_time_beg)/num_segs
  run_time_end = run_time_beg + run_time*segment
  run_time_beg = run_time_beg + run_time*(segment - 1)

  run_time_days_beg = int (run_time_beg/86400.0)
  run_time_secs_beg = nint(run_time_beg - run_time_days_beg*86400.0)
  
  Time_init = set_time ( run_time_secs_beg, run_time_days_beg )
  
  run_time_days_end = int (run_time_end/86400.0)
  run_time_secs_end = nint(run_time_end - run_time_days_end*86400.0)
         
  Time_end  = set_time ( run_time_secs_end, run_time_days_end )

  call get_time ( Time_end - Time_init, nsecs, ndays )
  run_time = real(nsecs) + 86400.0*ndays
  
  return
end subroutine get_times

!#######################################################################

subroutine write_timestamps ( Time_init, Time_end, num_segs )
    
type (time_type),                                 intent(inout) ::     &
                                                  Time_init, Time_end

integer, intent(in)                                             ::     &
                                                             num_segs

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

character(len=9) :: mname
character(len=4) :: yname='0000'
integer :: unit, nseg
integer :: year, month, day, hour, minute, second

!-----------------------------------------------------------------------
! write time stamps
!-----------------------------------------------------------------------

  unit = open_file ('time_stamp.out', action='append')

! overall beginning time
  call get_date ( Time_init, year, month, day, hour, minute, second )
  mname = month_name(month)
  if (mpp_pe() ==  mpp_root_pe()) write (unit,'(6i4,2x,a3)')           &
                     year, month, day, hour, minute, second, mname(1:3)
! overall ending time
  call get_date ( Time_end, year, month, day, hour, minute, second )
  mname = month_name(month)
  if (mpp_pe() ==  mpp_root_pe()) write (unit,'(6i4,2x,a3)')           &
                     year, month, day, hour, minute, second, mname(1:3)

  if ( num_segs > 1 ) then
     do nseg=1,num_segs
        call get_times ( Time_init, Time_end, seg=nseg, nsegs=num_segs )

! segment beginning time

  call get_date ( Time_init, year, month, day, hour, minute, second )
  mname = month_name(month)
  if (mpp_pe() ==  mpp_root_pe()) write (unit,'(6i4,2x,a3)')           &
                     year, month, day, hour, minute, second, mname(1:3)

! segment run time

  call get_time ( Time_end - Time_init, second, day )
  if (mpp_pe() ==  mpp_root_pe()) write (unit,'(2i10)') day, second
if (mpp_root_pe() == 0) print *,year, month, day, hour, minute, second, mname(1:3)
if (mpp_root_pe() == 0) print *,nseg,num_segs,day,second
     enddo
  endif

  call close_file (unit)

  return
end subroutine write_timestamps

!#######################################################################

subroutine error_handler ( message ) 
character(len=*), intent(in) :: message
   
  call error_mesg ('zetac_times_mod', message, FATAL)

  return
end subroutine error_handler

!#######################################################################

end module zetac_times_mod


