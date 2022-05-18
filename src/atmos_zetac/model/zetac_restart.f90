module zetac_restart_mod

!-----------------------------------------------------------------------
! module information
!-----------------------------------------------------------------------

use mpp_mod,                      only : mpp_pe, mpp_root_pe,          &
                                         stdout, stdlog,               &
                                         input_nml_file
use fms_mod,                      only : write_version_number,         &
                                         error_mesg, close_file,       &
                                         open_namelist_file,           &
                                         file_exist, check_nml_error,  &
                                         FATAL, WARNING
use mpp_io_mod,                   only : MPP_NETCDF, MPP_MULTI,        &
                                         MPP_SINGLE, MPP_RDONLY,       &
                                         MPP_OVERWR,                   &
                                         mpp_read, mpp_write,          &
                                         mpp_write_meta,               &
                                         mpp_open, mpp_close,          &
                                         mpp_get_axes,                 &
                                         mpp_get_fields,               &
                                         mpp_get_atts,                 & 
                                         mpp_get_axis_data,            &
                                         mpp_get_info, mpp_get_times,  &
                                         axistype, fieldtype
use diag_manager_mod,             only : get_base_date
use field_manager_mod,            only : model_atmos
use time_manager_mod,             only : time_type,                    &
                                         get_time, set_date, set_time, &
                                         operator(+), operator(-),     &
                                         operator(.lt.)
use get_cal_time_mod,             only : get_cal_time
use tracer_manager_mod,           only : get_tracer_names,             &
                                         get_tracer_index

use zetac_axis_names_mod,         only : naxes, xu_, xv_, yu_, yv_,    &
                                                xm_, ym_, zw_, zm_,    &
                                                yg_, t_
use zetac_axes_mod,               only : axes_restart
use zetac_field_names_mod
use zetac_horiz_grid_type_mod,    only : horiz_grid_type
use zetac_horiz_metric_type_mod,  only : horiz_metric_type
use zetac_ncdf_io_mod,            only : ncdf_read, ncdf_write
use zetac_time_pointers_mod,      only : ntime
use zetac_tracer_mod,             only : get_info
use zetac_update_halos_mod,       only : update_halos
use zetac_vert_metric_type_mod,   only : vert_metric_type

!-----------------------------------------------------------------------
! public interfaces
!-----------------------------------------------------------------------

implicit none

private
public read_restart, write_restart,                                    &
       read_restart_init, write_restart_init,                          &
       threading_rd, fileset_rd, threading_wt, fileset_wt

character(len=*), parameter :: module='zetac_restart_mod'

character*64, parameter :: filename_read = 'INPUT/zetac'
character*64, parameter :: filename_write = 'RESTART/zetac'
character*64, parameter :: filename_read_trac = 'INPUT/tracer.res.nc'
character*64, parameter :: filename_write_trac = 'RESTART/tracer.res'

integer  :: unit_read
integer  :: unit_write

integer :: nrefs=0
integer :: nprog=0, nprog2=0, ntend=0, nmean=0, ndiag=0, ndiag2=0

type(axistype) :: Axes_res(naxes)
type(axistype) :: Axes(naxes)

type(fieldtype), allocatable, dimension(:),  save :: Core_Fields
type(fieldtype), dimension(num_refs_fields), save :: Refs_Fields
type(fieldtype), dimension(num_stat_fields), save :: Stat_Fields
type(fieldtype), dimension(num_prog_fields), save :: Prog_Fields
type(fieldtype), dimension(num_prog_fields), save :: Prog_Fields2
type(fieldtype), dimension(num_tend_fields+2), save :: Tend_Fields
type(fieldtype), dimension(num_mean_fields+2), save :: Mean_Fields
type(fieldtype), dimension(num_diag_fields), save :: Diag_Fields
type(fieldtype), dimension(num_diag_fields), save :: Diag_Fields2

type(time_type), save :: Time_base
type(time_type), save :: Time_off

real, allocatable, dimension(:) :: xu
real, allocatable, dimension(:) :: yv
real, allocatable, dimension(:) :: zw

logical :: do_init=.true.
logical :: cold_read, cold_write=.true.

integer :: threading_rd, fileset_rd, threading_wt, fileset_wt 
integer :: secs, days, year_base, month_base, day_base, date(6)

!-----------------------------------------------------------------------
! version control information
!-----------------------------------------------------------------------

character(len=128)  :: version =  '$Id: zetac_restart.f90,v 1.1.2.8.2.10 2005/08/06 23:59:33 stg Exp $'
character(len=128)  :: tag     =  '$Name:  $'

contains

!#######################################################################

subroutine read_restart (                                              &
                                                                Mgrid, &
                                                           topog, sst, &
                                           ppref, uuref, taref, qvref, &
                                       ps, dp, uu, vv, oo, ta, qv, qc, &
                             utend, vtend, ttend, qtend, ttndf, qtndf, &
                               umean, vmean, tmean, qmean, udiv, vdiv, &
                                                        past, present )

!-----------------------------------------------------------------------
! reads time-dependent data from netcdf restart file
!-----------------------------------------------------------------------

type (horiz_grid_type),                           intent (in)     ::   &
                                                                Mgrid

real, dimension(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed),              &
                                                  intent (out)    ::   &
                                                           topog, sst

real, dimension(Mgrid%jbg:Mgrid%jeg, Mgrid%iz),   intent (out)    ::   &
                                           ppref, uuref, taref, qvref

real, dimension(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed,               &
                ntime),                           intent (inout)  ::   &
                                                                   ps

real, dimension(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed,               &
                Mgrid%kbd:Mgrid%ked, ntime),      intent (inout)  ::   &
                                           dp, uu, vv, oo, ta, qv, qc

real, dimension(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed,               &
                Mgrid%kbd:Mgrid%ked),             intent (out)    ::   &
                             utend, vtend, ttend, qtend, ttndf, qtndf

real, dimension(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed,               &
                Mgrid%kbd:Mgrid%ked),             intent (out)    ::   &
                               umean, vmean, tmean, qmean, udiv, vdiv

integer,                                          intent (in)     ::   &
                                                        past, present

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

character(len=16) :: name
integer :: i, ibd, ied
integer :: j, jbd, jed
integer :: k, kbd, ked
integer :: n, m, nt

!-----------------------------------------------------------------------
! read reference variable arrays with global axes
!-----------------------------------------------------------------------

  call mpp_read ( unit_read, Refs_fields(ppref_), ppref )
  call mpp_read ( unit_read, Refs_fields(uuref_), uuref )
  call mpp_read ( unit_read, Refs_fields(taref_), taref )
  call mpp_read ( unit_read, Refs_fields(qvref_), qvref )

  ibd = Mgrid%ibd
  ied = Mgrid%ied
  jbd = Mgrid%jbd
  jed = Mgrid%jed
  kbd = Mgrid%kbd
  ked = Mgrid%ked

!-----------------------------------------------------------------------
! read static variable arrays
!-----------------------------------------------------------------------

  do n=1,num_stat_fields
     call mpp_get_atts ( Stat_fields(n), name )
     if ( name == stat_names(topog_) ) then
        call ncdf_read ( Mgrid, unit_read, Stat_fields(n),             &
                                                   kbd, kbd, topog )
     endif
     call mpp_get_atts ( Stat_fields(n), name )
     if ( name == stat_names(sst_) ) then
        call ncdf_read ( Mgrid, unit_read, Stat_fields(n),             &
                                                     kbd, kbd, sst )
     endif
  enddo

  call update_halos ( Mgrid, topog )
  call update_halos ( Mgrid, sst )

!-----------------------------------------------------------------------
! read core variables at past time level
!-----------------------------------------------------------------------

  write (stdout(), '(/,"zetac_restart: reading variables" )')

  do n=1,num_prog_fields
     call mpp_get_atts ( Prog_fields(n), name )
     if ( name == prog_names(ps_) ) then
        call ncdf_read ( Mgrid, unit_read, Prog_fields(n),             &
                                               kbd, kbd, ps(:,:,past) )
     else if ( name == prog_names(uu_) ) then
        call ncdf_read ( Mgrid, unit_read, Prog_fields(n),             &
                                             kbd, ked, uu(:,:,:,past) )
     else if ( name == prog_names(vv_) ) then
        call ncdf_read ( Mgrid, unit_read, Prog_fields(n),             &
                                             kbd, ked, vv(:,:,:,past) )
     else if ( name == prog_names(oo_) ) then
        call ncdf_read ( Mgrid, unit_read, Prog_fields(n),             &
                                             kbd, ked, oo(:,:,:,past) )
     else if ( name == prog_names(ta_) ) then
        call ncdf_read ( Mgrid, unit_read, Prog_fields(n),             &
                                             kbd, ked, ta(:,:,:,past) )
     else if ( name == prog_names(qv_) ) then
        call ncdf_read ( Mgrid, unit_read, Prog_fields(n),             &
                                             kbd, ked, qv(:,:,:,past) )
     else if ( name == prog_names(qc_) ) then
        call ncdf_read ( Mgrid, unit_read, Prog_fields(n),             &
                                             kbd, ked, qc(:,:,:,past) )
     endif
  enddo
  
  do n=1,num_diag_fields
     call mpp_get_atts ( Diag_fields(n), name )
     if ( name == diag_names(dp_) ) then
        call ncdf_read ( Mgrid, unit_read, Diag_fields(n),             &
                                             kbd, ked, dp(:,:,:,past) )
     endif
  enddo

  call update_halos ( Mgrid, ps(:,:,past) ) 
  call update_halos ( Mgrid, uu(:,:,:,past) ) 
  call update_halos ( Mgrid, vv(:,:,:,past) )
  call update_halos ( Mgrid, oo(:,:,:,past) )
  call update_halos ( Mgrid, ta(:,:,:,past) )
  call update_halos ( Mgrid, qv(:,:,:,past) )
  call update_halos ( Mgrid, qc(:,:,:,past) )
  call update_halos ( Mgrid, dp(:,:,:,past) )
  ttndf = 0.
  qtndf = 0.
  udiv = 0.
  vdiv = 0.

  if ( cold_read ) then

     do j=jbd,jed
        do i=ibd,ied
          ps(i,j,present) = ps(i,j,past)
        enddo
     enddo
     do k=kbd,ked
        do j=jbd,jed
           do i=ibd,ied
              uu(i,j,k,present) = uu(i,j,k,past)
              vv(i,j,k,present) = vv(i,j,k,past)
              oo(i,j,k,present) = oo(i,j,k,past)
              ta(i,j,k,present) = ta(i,j,k,past)
              qv(i,j,k,present) = qv(i,j,k,past)
              qc(i,j,k,present) = qc(i,j,k,past)
              dp(i,j,k,present) = dp(i,j,k,past)
              utend(i,j,k) = 0.
              vtend(i,j,k) = 0.
              ttend(i,j,k) = 0.
              qtend(i,j,k) = 0.
              umean(i,j,k) = uu(i,j,k,present)
              vmean(i,j,k) = vv(i,j,k,present)
              tmean(i,j,k) = ta(i,j,k,present)
              qmean(i,j,k) = qv(i,j,k,present)
           enddo
        enddo
     enddo

  else

!-----------------------------------------------------------------------
! read core variables at present time level
!-----------------------------------------------------------------------

     write (stdout(), '("zetac_restart: reading time2 variables")' )

     do n=1,num_prog_fields
        call mpp_get_atts ( Prog_fields2(n), name )
        name = name(1:len(trim(name))-1)
        if ( name == prog_names(ps_) ) then
           call ncdf_read ( Mgrid, unit_read, Prog_fields2(n),         &
                                            kbd, kbd, ps(:,:,present) )
        else if ( name == prog_names(uu_) ) then
           call ncdf_read ( Mgrid, unit_read, Prog_fields2(n),         &
                                          kbd, ked, uu(:,:,:,present) )
        else if ( name == prog_names(vv_) ) then
           call ncdf_read ( Mgrid, unit_read, Prog_fields2(n),         &
                                          kbd, ked, vv(:,:,:,present) )
        else if ( name == prog_names(oo_) ) then
           call ncdf_read ( Mgrid, unit_read, Prog_fields2(n),         &
                                          kbd, ked, oo(:,:,:,present) )
        else if ( name == prog_names(ta_) ) then
           call ncdf_read ( Mgrid, unit_read, Prog_fields2(n),         &
                                          kbd, ked, ta(:,:,:,present) )
        else if ( name == prog_names(qv_) ) then
           call ncdf_read ( Mgrid, unit_read, Prog_fields2(n),         &
                                          kbd, ked, qv(:,:,:,present) )
        else if ( name == prog_names(qc_) ) then
           call ncdf_read ( Mgrid, unit_read, Prog_fields2(n),         &
                                          kbd, ked, qc(:,:,:,present) )
        endif
     enddo

     do n=1,num_diag_fields
        call mpp_get_atts ( Diag_fields2(n), name )
        name = name(1:len(trim(name))-1)
        if ( name == diag_names(dp_) ) then
           call ncdf_read ( Mgrid, unit_read, Diag_fields2(n),         &
                                          kbd, ked, dp(:,:,:,present) )
        endif
     enddo

!-----------------------------------------------------------------------
! read tendency fields
!-----------------------------------------------------------------------

!     do n=1,num_tend_fields+2
     do n=1,ntend
        if ( n == utend_ ) then
           call ncdf_read ( Mgrid, unit_read, Tend_fields(n),          &
                                               kbd, ked, utend(:,:,:) )
        else if ( n == vtend_ ) then
           call ncdf_read ( Mgrid, unit_read, Tend_fields(n),          &
                                               kbd, ked, vtend(:,:,:) )
        else if ( n == ttend_ ) then
           call ncdf_read ( Mgrid, unit_read, Tend_fields(n),          &
                                               kbd, ked, ttend(:,:,:) )
        else if ( n == qtend_ ) then
           call ncdf_read ( Mgrid, unit_read, Tend_fields(n),          &
                                               kbd, ked, qtend(:,:,:) )
        else if (n == num_tend_fields+1) then
           call ncdf_read ( Mgrid, unit_read, Tend_fields(n),          &
                                               kbd, ked, ttndf(:,:,:) )
        else if (n == num_tend_fields+2) then
           call ncdf_read ( Mgrid, unit_read, Tend_fields(n),          &
                                               kbd, ked, qtndf(:,:,:) )
        endif
     enddo

!-----------------------------------------------------------------------
! read mean fields
!-----------------------------------------------------------------------

!     do n=1,num_mean_fields+2
     do n=1,nmean
        if ( n == umean_ ) then
           call ncdf_read ( Mgrid, unit_read, Mean_fields(n),          &
                                               kbd, ked, umean(:,:,:) )
        else if ( n == vmean_ ) then
           call ncdf_read ( Mgrid, unit_read, Mean_fields(n),          &
                                               kbd, ked, vmean(:,:,:) )
        else if ( n == tmean_ ) then
           call ncdf_read ( Mgrid, unit_read, Mean_fields(n),          &
                                               kbd, ked, tmean(:,:,:) )
        else if ( n == qmean_ ) then
           call ncdf_read ( Mgrid, unit_read, Mean_fields(n),          &
                                               kbd, ked, qmean(:,:,:) )
        else if ( n == num_mean_fields+1 ) then
           call ncdf_read ( Mgrid, unit_read, Mean_fields(n),          &
                                                kbd, ked, udiv(:,:,:) )
        else if ( n == num_mean_fields+2 ) then
           call ncdf_read ( Mgrid, unit_read, Mean_fields(n),          &
                                                kbd, ked, vdiv(:,:,:) )
        endif
     enddo

     call update_halos ( Mgrid, ps(:,:,present) ) 
     call update_halos ( Mgrid, uu(:,:,:,present) ) 
     call update_halos ( Mgrid, vv(:,:,:,present) ) 
     call update_halos ( Mgrid, oo(:,:,:,present) ) 
     call update_halos ( Mgrid, ta(:,:,:,present) )
     call update_halos ( Mgrid, qv(:,:,:,present) )
     call update_halos ( Mgrid, qc(:,:,:,present) )
     call update_halos ( Mgrid, dp(:,:,:,present) )

     call update_halos ( Mgrid, ttend(:,:,:) )
     call update_halos ( Mgrid, qtend(:,:,:) )
     call update_halos ( Mgrid, umean(:,:,:) )
     call update_halos ( Mgrid, vmean(:,:,:) )
     call update_halos ( Mgrid, tmean(:,:,:) )
     call update_halos ( Mgrid, qmean(:,:,:) )
     call update_halos ( Mgrid, udiv(:,:,:) )
     call update_halos ( Mgrid, vdiv(:,:,:) )

  endif

!-----------------------------------------------------------------------
! close input restart file(s)
!-----------------------------------------------------------------------

!!$  call mpp_close ( unit_read_core )                                  !stg
!!$  if ( unit_read_save > 0 ) call mpp_close ( unit_read_save )

  return
end subroutine read_restart

!#######################################################################

subroutine write_restart (                                             &
                                               Mgrid, Time, Time_step, &
                                                           topog, sst, &
                                           ppref, uuref, taref, qvref, &
                                       ps, dp, uu, vv, oo, ta, qv, qc, &
                             utend, vtend, ttend, qtend, ttndf, qtndf, &
                               umean, vmean, tmean, qmean, udiv, vdiv, &
                                                       ltime1, ltime2 )

!-----------------------------------------------------------------------
! writes time-dependent data to netcdf restart file
!-----------------------------------------------------------------------

type (horiz_grid_type)                                         ::      &
                                                                Mgrid

type (time_type)                                               ::      &
                                                      Time, Time_step

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed),             &
                                                intent (in)    ::      &
                                                           topog, sst

real, dimension(Mgrid%jbg:Mgrid%jeg, Mgrid%iz), intent (in)    ::      &
                                           ppref, uuref, taref, qvref

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                ntime),                         intent (in)    ::      &
                                                                   ps

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked, ntime),    intent (in)    ::      &
                                           dp, uu, vv, oo, ta, qv, qc

real, optional, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,    &
                Mgrid%kbd:Mgrid%ked),           intent (in)    ::      &
                             utend, vtend, ttend, qtend, ttndf, qtndf

real, optional, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,    &
                Mgrid%kbd:Mgrid%ked),           intent (in)    ::      &
                               umean, vmean, tmean, qmean, udiv, vdiv

integer, optional,                              intent (in)    ::      &
                                                       ltime1, ltime2

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

integer :: kbd, ked
integer :: n, ltime

real    :: time_now, time_end

  if ( Time .lt. Time_base )                                           &
      call error_handler ('restart time earlier than base time', FATAL)

  kbd = Mgrid%kbd
  ked = Mgrid%ked

!-----------------------------------------------------------------------
! write reference variables with global axes
!-----------------------------------------------------------------------

  if ( mpp_pe() == mpp_root_pe() ) then
     call mpp_write ( unit_write, Refs_fields(ppref_), ppref )  
     call mpp_write ( unit_write, Refs_fields(uuref_), uuref )  
     call mpp_write ( unit_write, Refs_fields(taref_), taref )
     call mpp_write ( unit_write, Refs_fields(qvref_), qvref )
  endif

!-----------------------------------------------------------------------
! write static variables
!-----------------------------------------------------------------------

  call ncdf_write ( Mgrid, unit_write, Stat_fields(topog_),            &
                                                      kbd, kbd, topog )
  call ncdf_write ( Mgrid, unit_write, Stat_fields(sst_),              &
                                                        kbd, kbd, sst )

  call get_time (Time - Time_base, secs, days)
  time_end = real(secs)/60.0 + 1440.0*days

  if ( cold_write ) then
     ltime = 1
     time_now = time_end
  else
     ltime = ltime1
     call get_time (Time_step, secs)
     time_now = time_end - real(secs)/60.0
  endif

!-----------------------------------------------------------------------
! write past time level
!-----------------------------------------------------------------------

  write (stdout(), '(/,"zetac_restart: writing time = ",f13.2)')       &
                                                             time_now

  do n=1,num_prog_fields
     if ( n == ps_ ) then
        call ncdf_write ( Mgrid, unit_write, Prog_fields(n),           &
                                    kbd, kbd, ps(:,:,ltime) )
     else if ( n == uu_ ) then
        call ncdf_write ( Mgrid, unit_write, Prog_fields(n),           &
                                  kbd, ked, uu(:,:,:,ltime) )
     else if ( n == vv_ ) then
        call ncdf_write ( Mgrid, unit_write, Prog_fields(n),           &
                                  kbd, ked, vv(:,:,:,ltime) )
     else if ( n == oo_ ) then
        call ncdf_write ( Mgrid, unit_write, Prog_fields(n),           &
                                  kbd, ked, oo(:,:,:,ltime) )
     else if ( n == ta_ ) then
        call ncdf_write ( Mgrid, unit_write, Prog_fields(n),           &
                                  kbd, ked, ta(:,:,:,ltime) )
     else if ( n == qv_ ) then
        call ncdf_write ( Mgrid, unit_write, Prog_fields(n),           &
                                  kbd, ked, qv(:,:,:,ltime) )
     else if ( n == qc_ ) then
        call ncdf_write ( Mgrid, unit_write, Prog_fields(n),           &
                                  kbd, ked, qc(:,:,:,ltime) )
     endif
  enddo
  
  do n=1,num_diag_fields
     if ( n == dp_ ) then
        call ncdf_write ( Mgrid, unit_write, Diag_fields(n),           &
             kbd, ked, dp(:,:,:,ltime) )
     endif
  enddo

  if ( .not. cold_write ) then

!-----------------------------------------------------------------------
! write present time level
!-----------------------------------------------------------------------

     write (stdout(), '("zetac_restart: writing time = ",f13.2)')      &
                                                             time_end

     do n=1,num_prog_fields
        if ( n == ps_ ) then
           call ncdf_write ( Mgrid, unit_write, Prog_fields2(n),       &
                                   kbd, kbd, ps(:,:,ltime2) )
        else if ( n == uu_ ) then
           call ncdf_write ( Mgrid, unit_write, Prog_fields2(n),       &
                                 kbd, ked, uu(:,:,:,ltime2) )
        else if ( n == vv_ ) then
           call ncdf_write ( Mgrid, unit_write, Prog_fields2(n),       &
                                 kbd, ked, vv(:,:,:,ltime2) )
        else if ( n == oo_ ) then
           call ncdf_write ( Mgrid, unit_write, Prog_fields2(n),       &
                                 kbd, ked, oo(:,:,:,ltime2) )
        else if ( n == ta_ ) then
           call ncdf_write ( Mgrid, unit_write, Prog_fields2(n),       &
                                 kbd, ked, ta(:,:,:,ltime2) )
        else if ( n == qv_ ) then
           call ncdf_write ( Mgrid, unit_write, Prog_fields2(n),       &
                                 kbd, ked, qv(:,:,:,ltime2) )
        else if ( n == qc_ ) then
           call ncdf_write ( Mgrid, unit_write, Prog_fields2(n),       &
                                 kbd, ked, qc(:,:,:,ltime2) )
        endif
     enddo
     
     do n=1,num_diag_fields
        if ( n == dp_ ) then
           call ncdf_write ( Mgrid, unit_write, Diag_fields2(n),       &
                                 kbd, ked, dp(:,:,:,ltime2) )
        endif
     enddo

!-----------------------------------------------------------------------
! write tendency fields
!-----------------------------------------------------------------------

     do n=1,num_tend_fields+2
        if ( n == utend_ ) then
           call ncdf_write ( Mgrid, unit_write, Tend_fields(n),        &
                                               kbd, ked, utend(:,:,:) )
        else if ( n == vtend_ ) then
           call ncdf_write ( Mgrid, unit_write, Tend_fields(n),        &
                                               kbd, ked, vtend(:,:,:) )
        else if ( n == ttend_ ) then
           call ncdf_write ( Mgrid, unit_write, Tend_fields(n),        &
                                               kbd, ked, ttend(:,:,:) )
        else if ( n == qtend_ ) then
           call ncdf_write ( Mgrid, unit_write, Tend_fields(n),        &
                                               kbd, ked, qtend(:,:,:) )
        else if ( n == num_tend_fields+1 ) then
           call ncdf_write ( Mgrid, unit_write, Tend_fields(n),        &
                                               kbd, ked, ttndf(:,:,:) )
        else if ( n == num_tend_fields+2 ) then
           call ncdf_write ( Mgrid, unit_write, Tend_fields(n),        &
                                               kbd, ked, qtndf(:,:,:) )
        endif
     enddo

!-----------------------------------------------------------------------
! write mean fields
!-----------------------------------------------------------------------

     do n=1,num_mean_fields+2
        if ( n == umean_ ) then
           call ncdf_write ( Mgrid, unit_write, Mean_fields(n),        &
                                               kbd, ked, umean(:,:,:) )
        else if ( n == vmean_ ) then
           call ncdf_write ( Mgrid, unit_write, Mean_fields(n),        &
                                               kbd, ked, vmean(:,:,:) )
        else if ( n == tmean_ ) then
           call ncdf_write ( Mgrid, unit_write, Mean_fields(n),        &
                                               kbd, ked, tmean(:,:,:) )
        else if ( n == qmean_ ) then
           call ncdf_write ( Mgrid, unit_write, Mean_fields(n),        &
                                               kbd, ked, qmean(:,:,:) )
        else if ( n == num_mean_fields+1 ) then
           call ncdf_write ( Mgrid, unit_write, Mean_fields(n),        &
                                                kbd, ked, udiv(:,:,:) )
        else if ( n == num_mean_fields+2 ) then
           call ncdf_write ( Mgrid, unit_write, Mean_fields(n),        &
                                                kbd, ked, vdiv(:,:,:) )
        endif
     enddo

  endif  ! cold_write

!-----------------------------------------------------------------------
! close output restart file(s)
!-----------------------------------------------------------------------

  call mpp_close ( unit_write )

  return
end subroutine write_restart

!#######################################################################

subroutine read_restart_init (                                         &
                                               Mgrid, Time, Time_step, &
                                   rlonmin, rlonmax, rlatmin, rlatmax, &
                                                                 zeta, &
                                                            coldstart )


type (horiz_grid_type), intent (in)                           ::       &
                                                                Mgrid

type (time_type),       intent (in)                           ::       &
                                                      Time, Time_step

real,                   intent (in)                           ::       &
                                   rlonmin, rlonmax, rlatmin, rlatmax

real, dimension(Mgrid%iz),                                             &
                         intent (in)                         ::        &
                                                                 zeta

logical,                intent (out)                          ::       &
                                                            coldstart

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

character*64 :: units, calendar_type
character*64 :: file, msg
character*64 :: name, name2

real    :: dlon, dlon_rst, dlat, dlat_rst, time_rst, dt_rst
real    :: zw_rst(Mgrid%iz)

real    :: time_cal, time_bas, dt

integer :: ibg, ieg
integer :: jbg, jeg
integer :: ix, ix_rst=0
integer :: iy, iy_rst=0
integer :: iz, iz_rst=0
integer :: nt, nvar, ndim, natt
integer :: n, m, length

  ibg = Mgrid%ibg
  ieg = Mgrid%ieg
  jbg = Mgrid%jbg
  jeg = Mgrid%jeg
  ix  = Mgrid%ix
  iy  = Mgrid%iy
  iz  = Mgrid%iz

!-----------------------------------------------------------------------
! initialize module
!-----------------------------------------------------------------------

  if (do_init) then
     call restart_init ( Mgrid )
  endif

!-----------------------------------------------------------------------
! open input restart file and save fieldtype metadata
!-----------------------------------------------------------------------

  file = trim(filename_read) // '.res.nc'

  call mpp_open (                unit_read, file(1:len(trim(file))-3), &
                                   action=MPP_RDONLY, form=MPP_NETCDF, &
                           threading=threading_rd, fileset=fileset_rd )

  call mpp_get_info ( unit_read, ndim, nvar, natt, nt )
  allocate ( Core_Fields(nvar) )
  call mpp_get_fields ( unit_read, Core_Fields )

  do n=1,nvar
     call mpp_get_atts ( Core_fields(n), name )
     do m=1,num_stat_fields
        if ( name == stat_names(m) ) then
           nrefs = nrefs+1 ; Stat_fields(m) = Core_fields(n)
        endif
     enddo
     do m=1,num_refs_fields
        if ( name == refs_names(m) ) then
           nrefs = nrefs+1 ; Refs_fields(m) = Core_fields(n)
        endif
     enddo
     do m=1,num_prog_fields
        if ( name == prog_names(m) ) then
           nprog = nprog+1 ; Prog_fields(m) = Core_fields(n)
        endif
     enddo
     do m=1,num_prog_fields
        name2 = name(1:len(trim(name))-1)
        if ( name2 == prog_names(m) ) then
           nprog2 = nprog2+1 ; Prog_fields2(m) = Core_fields(n)
        endif
     enddo
     do m=1,num_tend_fields
        if ( name == tend_names(m) ) then
           ntend = ntend+1 ; Tend_fields(m) = Core_fields(n)
        endif
     enddo
     if ( name == 'ttndf' ) then
           ntend = ntend+1 ; Tend_fields(ntend) = Core_fields(n)
     else if ( name == 'qtndf' ) then
           ntend = ntend+1 ; Tend_fields(ntend) = Core_fields(n)
     endif
     do m=1,num_mean_fields
        if ( name == mean_names(m) ) then
           nmean = nmean+1 ; Mean_fields(m) = Core_fields(n)
        endif
     enddo
     if ( name == 'udiv' ) then
           nmean = nmean+1 ; Mean_fields(nmean) = Core_fields(n)
     else if ( name == 'vdiv' ) then
           nmean = nmean+1 ; Mean_fields(nmean) = Core_fields(n)
     endif
     do m=1,num_diag_fields
        if ( name == diag_names(m) ) then
           ndiag = ndiag+1 ; Diag_fields(m) = Core_fields(n)
        endif
     enddo
     do m=1,num_diag_fields
        name2 = name(1:len(trim(name))-1)
        if ( name2 == diag_names(m) ) then
           ndiag2 = ndiag2+1 ; Diag_fields2(m) = Core_fields(n)
        endif
     enddo
  enddo

  if ( nprog < num_prog_fields )          &
         call error_handler ('fields missing from restart file', FATAL)

  cold_read = ( nprog2 == 0 )

!-----------------------------------------------------------------------
! read axes and check for agreement with namelist
!-----------------------------------------------------------------------

  if ( ndim /= naxes-1 )                                               &
        call error_handler ('unexpected number of axes found', WARNING)

  call mpp_get_axes ( unit_read, Axes(1:ndim) )

  do n=1,ndim
     call mpp_get_atts ( Axes(n), name=name, units=units,              &
                                   len=length, calendar=calendar_type )
     if ( name == 'lon_u' ) then
        ix_rst = length - 4
        if ( ix_rst /= ix ) then
           write (msg,'("grid dimension ix/new,old/ =",2i6)') ix, ix_rst
           call error_handler ( msg, FATAL )
        endif
        allocate ( xu(ibg:ieg) )
        call mpp_get_axis_data ( Axes(n), xu )
     else if ( name == 'lat_v' ) then
        iy_rst = length - 4
        if ( iy_rst /= iy ) then
           write (msg,'("grid dimension iy/new,old/ =",2i6)') iy, iy_rst
           call error_handler ( msg, FATAL )
        endif
        allocate ( yv(jbg:jeg) )
        call mpp_get_axis_data ( Axes(n), yv )
     else if ( name == 'z_w' ) then
        iz_rst = length
        if ( iz_rst /= iz ) then
           write (msg,'("grid dimension iz/new,old/ =",2i6)') iz, iz_rst
           call error_handler ( msg, FATAL )
        endif
        allocate ( zw(iz_rst) )
        call mpp_get_axis_data ( Axes(n), zw )
     else if ( name == 'time' ) then
        time_cal = 0.0
        if ( len(trim(units)) > 7 .and. year_base /= 0 ) then
           Time_off = get_cal_time ( time_cal, units,                  &
                                                  trim(calendar_type) )
           call get_time (Time_off, secs, days)
           time_cal = real(secs)/60.0 + 1440.0*days
        endif
     endif
  enddo

!-----------------------------------------------------------------------
! check grid for agreement with namelist
!-----------------------------------------------------------------------

  if ( ix_rst == 0 .or. iy_rst == 0 .or. iz_rst == 0 )                 &
                  call error_handler ('required axis not found', FATAL)

  dlon = (rlonmax - rlonmin)/float(ix)
  dlon_rst = xu(2) - xu(1)
  if (abs(rlonmax - xu(ix)) > 0.01*dlon) then
     write ( msg, '("rlonmax/new,old/ =",2f9.3)' ) rlonmax, xu(ix)
     call error_handler ( msg, WARNING )
  endif
  if (abs(dlon - dlon_rst) > 0.01*dlon) then
     write ( msg, '("dlon/new,old/ =",2f9.3)' ) dlon, dlon_rst
     call error_handler ( msg, WARNING )
  endif

  dlat = (rlatmax - rlatmin)/float(iy)
  dlat_rst = yv(2) - yv(1)
  if (abs(rlatmax - yv(iy)) > 0.01*dlat) then
     write ( msg, '("rlatmax/new,old/ =", 2f9.3)' ) rlatmax, yv(iy)
     call error_handler ( msg, WARNING )
  endif
  if (abs(dlat - dlat_rst) > 0.01*dlat) then
     write ( msg, '("dlat/new,old/ =", 2f9.3)' ) dlat, dlat_rst
     call error_handler ( msg, WARNING )
  endif

  zw_rst = zw - zeta
  if (any(abs(zw_rst) > 0.0001)) then
     write ( msg, '("zeta/new,old/ =", 2f9.2)' ) zw(2), zeta(2)
     call error_handler ( msg, WARNING )
  endif

  cold_write = .false.
  coldstart = cold_read

  return
end subroutine read_restart_init

!#######################################################################

subroutine write_restart_init ( Mgrid, Hmetric, Vmetric )

type (horiz_grid_type),                           intent (in)    ::    &
                                                                Mgrid

type (horiz_metric_type),                         intent (in)    ::    &
                                                              Hmetric

type (vert_metric_type),                          intent (in)    ::    &
                                                              Vmetric

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

character*64 :: longname, units
character*64 :: file
character(len=16) :: name
integer :: m, n

type(axistype), dimension(3) :: Axes_uum, Axes_vvm, Axes_mmm, Axes_mmw
type(axistype), dimension(2) :: Axes_mm
type(axistype), dimension(2) :: Axes_mmg

!-----------------------------------------------------------------------
! initialize module
!-----------------------------------------------------------------------

  if (do_init) then
     call restart_init ( Mgrid )
  endif

!-----------------------------------------------------------------------
! open output restart files and register axes
!-----------------------------------------------------------------------

  file = trim(filename_write) // '.res'
  call mpp_open     (                                unit_write, file, &
                                   action=MPP_OVERWR, form=MPP_NETCDF, &
                           threading=threading_wt, fileset=fileset_wt )  

  call axes_restart (                         Mgrid, Hmetric, Vmetric, &
                                                 unit_write, Axes_res )!, &
!                       year=year_base, month=month_base, day=day_base )

!-----------------------------------------------------------------------
! define axis sets
!-----------------------------------------------------------------------

  Axes_uum = (/Axes_res(xu_),Axes_res(yu_),Axes_res(zm_)/)
  Axes_vvm = (/Axes_res(xv_),Axes_res(yv_),Axes_res(zm_)/)
  Axes_mmm = (/Axes_res(xm_),Axes_res(ym_),Axes_res(zm_)/)
  Axes_mmw = (/Axes_res(xm_),Axes_res(ym_),Axes_res(zw_)/)
  Axes_mm  = (/Axes_res(xm_),Axes_res(ym_)/)
  Axes_mmg = (/Axes_res(yg_),Axes_res(zm_)/)

!-----------------------------------------------------------------------
! register reference fields
!-----------------------------------------------------------------------

  do n=1,num_refs_fields
     call get_info ( refs_names(n), longname, units )
     call mpp_write_meta( unit_write, Refs_fields(n),                  & 
                     Axes_mmg, refs_names(n), units, longname, pack=1 )
  enddo

!-----------------------------------------------------------------------
! register static fields
!-----------------------------------------------------------------------

  do n=1,num_stat_fields
     call get_info ( stat_names(n), longname, units )
     call mpp_write_meta( unit_write, Stat_fields(n),                  &
                      Axes_mm, stat_names(n), units, longname, pack=1 )
  enddo

!-----------------------------------------------------------------------
! register prognostic fields
!-----------------------------------------------------------------------

  do n=1,num_prog_fields
     call get_info ( prog_names(n), longname, units )
     if ( n == ps_ ) then
        call mpp_write_meta( unit_write, Prog_fields(n),               &
               Axes_mm, prog_names(n), units, longname, pack=1 )
     else if ( n == uu_ ) then
        call mpp_write_meta( unit_write, Prog_fields(n),               &
               Axes_uum, prog_names(n), units, longname, pack=1 )
     else if ( n == vv_ ) then
        call mpp_write_meta( unit_write, Prog_fields(n),               &
               Axes_vvm, prog_names(n), units, longname, pack=1 )
     else if ( n == oo_ ) then
        call mpp_write_meta( unit_write, Prog_fields(n),               &
               Axes_mmw, prog_names(n), units, longname, pack=1 )
     else if ( n == ta_ .or. n == qv_ .or. n == qc_ ) then
        call mpp_write_meta( unit_write, Prog_fields(n),               &
               Axes_mmm, prog_names(n), units, longname, pack=1 )
     endif
  enddo

  do n=1,num_diag_fields
     call get_info ( diag_names(n), longname, units )
     if ( n == dp_ ) then
        call mpp_write_meta( unit_write, Diag_fields(n),               &
               Axes_mmm, diag_names(n), units, longname, pack=1 )
     endif
  enddo

!-----------------------------------------------------------------------
! register present-time and tendency fields
!-----------------------------------------------------------------------

  if (.not. cold_write) then

     do n=1,num_prog_fields
        call get_info ( prog_names(n), longname, units )
        name = trim(prog_names(n))//"2"
        if ( n == ps_ ) then
           call mpp_write_meta( unit_write, Prog_fields2(n),           &
                  Axes_mm, name, units, longname, pack=1 )
        else if ( n == uu_ ) then
           call mpp_write_meta( unit_write, Prog_fields2(n),           &
                  Axes_uum, name, units, longname, pack=1 )
        else if ( n == vv_ ) then
           call mpp_write_meta( unit_write, Prog_fields2(n),           &
                  Axes_vvm, name, units, longname, pack=1 )
        else if ( n == oo_ ) then
           call mpp_write_meta( unit_write, Prog_fields2(n),           &
                  Axes_mmw, name, units, longname, pack=1 )
        else if ( n == ta_ .or. n == qv_ .or. n == qc_ ) then
           call mpp_write_meta( unit_write, Prog_fields2(n),           &
                  Axes_mmm, name, units, longname, pack=1 )
        endif
     enddo

! present-time fields
     
     do n=1,num_diag_fields
        call get_info ( diag_names(n), longname, units )
        name = trim(diag_names(n))//"2"
        if ( n == dp_ ) then
           call mpp_write_meta( unit_write, Diag_fields2(n),               &
               Axes_mmm, name, units, longname, pack=1 )
        endif
     enddo

! tendency fields

     do n=1,num_tend_fields
        if (tend_names(n) == 'utend') then
!           longname = 'zonal velocity tendency'
           longname = 'zonal velocity tendency from time filter'
           units = 'm/s2'
        else if (tend_names(n) == 'vtend') then
!           longname = 'merid velocity tendency'
           longname = 'merid velocity tendency from time filter'
           units = 'm/s2'
        else if (tend_names(n) == 'ttend') then
           longname = 'temperature tendency'
           units = 'K/s'
        else if (tend_names(n) == 'qtend') then
           longname = 'specific humidity tendency'
           units = 'kg/kg/s'
        endif
        call mpp_write_meta( unit_write, Tend_fields(n),               &
               Axes_mmm, tend_names(n), units, longname, pack=1 )
     enddo
     n = num_tend_fields+1
     name = 'ttndf'
     longname = 'temperature tendency from time filter'
     units = 'K/s'
     call mpp_write_meta( unit_write, Tend_fields(n),                  &
                        Axes_mmm, name, units, longname, pack=1 )
     n = n+1
     name = 'qtndf'
     longname = 'specific humidity tendency from time filter'
     units = 'kg/kg/s'
     call mpp_write_meta( unit_write, Tend_fields(n),                  &
                        Axes_mmm, name, units, longname, pack=1 )

! mean fields

     do n=1,num_mean_fields
        if (mean_names(n) == 'umean') then
           longname = 'smooth zonal veloctity'
           units = 'm/s'
           call mpp_write_meta( unit_write, Mean_fields(n),            &
               Axes_uum, mean_names(n), units, longname, pack=1 )
        else if (mean_names(n) == 'vmean') then
           longname = 'smooth meridional veloctity'
           units = 'm/s'
           call mpp_write_meta( unit_write, Mean_fields(n),            &
               Axes_vvm, mean_names(n), units, longname, pack=1 )
        else if (mean_names(n) == 'tmean') then
           longname = 'smooth temperature'
           units = 'K'
           call mpp_write_meta( unit_write, Mean_fields(n),            &
               Axes_mmm, mean_names(n), units, longname, pack=1 )
        else if (mean_names(n) == 'qmean') then
           longname = 'smooth specific humidity'
           units = 'kg/kg'
           call mpp_write_meta( unit_write, Mean_fields(n),            &
               Axes_mmm, mean_names(n), units, longname, pack=1 )
        endif
     enddo
     n = num_mean_fields+1
     name = 'udiv'
     longname = 'divergent zonal velocity'
     units = 'm/s'
     call mpp_write_meta( unit_write, Mean_fields(n),                  &
                        Axes_uum, name, units, longname, pack=1 )
     n = n+1
     name = 'vdiv'
     longname = 'divergent meridional velocity'
     units = 'm/s'
     call mpp_write_meta( unit_write, Mean_fields(n),                  &
                        Axes_vvm, name, units, longname, pack=1 )

  endif  ! test cold_write

!-----------------------------------------------------------------------
! write axis info (tracer axes are written in atmosphere_mod)
!-----------------------------------------------------------------------

  do m=1,naxes-1
     call mpp_write ( unit_write, Axes_res(m) )
  enddo

  return
end subroutine write_restart_init

!#######################################################################

subroutine restart_init ( Mgrid )

type (horiz_grid_type),                           intent (in)    ::    &
                                                                Mgrid

integer :: io, unit_nml, ierr

character(len=6) :: threading_read='multi', fileset_read='single',      &
                    threading_write='single', fileset_write='single'

namelist /fms_io_nml/ threading_read

!-----------------------------------------------------------------------
! read namelist
!-----------------------------------------------------------------------

  read (input_nml_file, nml=fms_io_nml, iostat=io)
  ierr = check_nml_error(io,'fms_io_nml')

!-----------------------------------------------------------------------
! write version number and namelist to logfile
!-----------------------------------------------------------------------

  call write_version_number (version, tag)
  if (mpp_pe() == mpp_root_pe())                                       &
                                       write (stdlog(), nml=fms_io_nml)

  if ( threading_read  == 'single' ) then
     threading_rd = MPP_SINGLE
  else
     threading_rd = MPP_MULTI
  endif
  if ( fileset_read  == 'single' ) then
     fileset_rd = MPP_SINGLE
  else
     fileset_rd = MPP_MULTI
  endif
  if ( threading_write  == 'single' ) then
     threading_wt = MPP_SINGLE
  else
     threading_wt = MPP_MULTI
  endif
  if ( fileset_write  == 'single' ) then
     fileset_wt = MPP_SINGLE
  else
     fileset_wt = MPP_MULTI
  endif

  call get_base_date ( date(1), date(2), date(3),                      &
                       date(4), date(5), date(6) )

  year_base = date(1)
  month_base = date(2)
  day_base = date(3)

  if ( year_base == 0 ) then
     Time_base = set_time ( (date(4)*60 + date(5))*60 + date(6),       &
                             date(3) )
  else
     Time_base = set_date ( date(1), date(2), date(3),                 &
                            date(4), date(5), date(6) )
  endif
  
  do_init = .false.

  return
end subroutine restart_init

!#######################################################################

subroutine error_handler ( message, level ) 
character(len=*), intent(in) :: message
integer,          intent(in) :: level

  if ( level == FATAL .or. mpp_pe() == mpp_root_pe() ) then
     call error_mesg (module, message, level)
  endif

  return
end subroutine error_handler

!#######################################################################

end module zetac_restart_mod

