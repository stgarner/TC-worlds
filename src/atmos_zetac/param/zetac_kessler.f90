module zetac_kessler_mod

!-----------------------------------------------------------------------
! module information
!-----------------------------------------------------------------------

use mpp_mod,            only : mpp_pe, mpp_root_pe, input_nml_file
use mpp_io_mod,         only : axistype, fieldtype
use fms_mod,            only : file_exist, check_nml_error,            &
                               error_mesg, FATAL,                      &
                               close_file, open_namelist_file,         &
                               stdlog, write_version_number
use constants_mod,      only : cp_air, hlv, hlf, tfreeze, grav

use diag_manager_mod,            only : register_diag_field
use field_manager_mod,           only : model_atmos
use time_manager_mod,            only : time_type
use tracer_manager_mod,          only : get_tracer_names,              &
                                        get_tracer_index,              &
                                        NO_TRACER

use zetac_axes_mod,              only : axid, gridm
use zetac_axis_names_mod,        only : xm_, ym_, zm_
use zetac_convert_var_mod,       only : get_pp
use zetac_extrap_var_mod,        only : extrap_tracer
use zetac_field_names_mod,       only : qc_
use zetac_horiz_grid_type_mod,   only : horiz_grid_type
use zetac_horiz_metric_type_mod, only : horiz_metric_type
use zetac_moisture_mod,          only : get_qsat, get_dlnesatdt
use zetac_ncdf_io_mod,           only : ncdf_fms_write
use zetac_phys_con_mod,          only : tnuke, eps ! Rd/Rv
use zetac_update_halos_mod,      only : update_halos
use zetac_vert_metric_type_mod,  only : vert_metric_type

!-----------------------------------------------------------------------
! public interfaces
!-----------------------------------------------------------------------

implicit none
private
public kessler_mp, kessler_init, kessler_end

character (len=*),  parameter :: module='zetac_kessler_mod'
character (len=16), parameter :: model='atmos_mod'
character (len=16) :: name, units
character (len=64) :: longname

real, parameter :: tsleet=tfreeze+2.0
real, parameter :: silly=-9999.9
real, parameter :: tiny=1.e-10

real, parameter :: eaccr=0.875   ! exponent for accretion
real, parameter :: eevap=0.525   ! exponent for evaporation

integer :: id_rain_mp, id_snow_mp, id_tdt_mp, id_qdt_mp, id_ldt_mp

real, allocatable, dimension(:,:,:) :: hlat, auto, wice, qsat, dqdt
real, allocatable, dimension(:,:,:) :: d_ta, d_qv, d_qc
real, allocatable, dimension(:,:,:) :: tabs1, qvap1, qcld1
             
! namelist:
real :: revap=0.2e-2, raccr=1.e-2
real :: rauto=0.2e-2, cthresh=5.e-4, rthresh=5.e-4
logical :: do_evap=.false., do_auto=.false., do_ice=.false.

integer :: ibc, iec
integer :: jbc, jec
integer :: kbd, ked

!-----------------------------------------------------------------------
! version control information
!-----------------------------------------------------------------------

character(len=128)  :: version =  '$Id: zetac_kessler.f90,v 1.1.2.9.2.4.2.5 2005/07/18 18:44:25 stg Exp $'
character(len=128)  :: tag     =  '$Name:  $'

contains

!#######################################################################

subroutine kessler_mp (                                                &
                                                          Mgrid, Time, &
                                                            pres, dpm, &
                                                     tabs, qvap, qcld, &
                                                           rain, snow, &
                                                                 delt )

type (horiz_grid_type),                      intent (in)    ::         &
                                                                Mgrid
       
type (time_type),                            intent (in)    ::         &
                                                                 Time

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),        intent (in)    ::         &
                                                            pres, dpm

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),        intent (inout) ::         &
                                                     tabs, qvap, qcld

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed),             &
                                             intent (inout)   ::         &
                                                           rain, snow

real,                                        intent (in)    ::         &
                                                                 delt

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

real    :: dt, wsleet
integer :: i, j, k

!-----------------------------------------------------------------------
! factor for diagnostic ice phase
!-----------------------------------------------------------------------

  if ( do_ice ) then
     wice = tanh(dim(tfreeze - tnuke, tabs)/tnuke)
     hlat = hlv + hlf*wice
  else
     wice = 0.0
     hlat = hlv
  endif

!-----------------------------------------------------------------------
! get changes in qv, qc and ta due to condensation
!-----------------------------------------------------------------------

  tabs1 = tabs
  qvap1 = qvap
  qcld1 = qcld

  call get_qsat ( pres, tabs, qsat, wice )

  call condense ( Mgrid, pres, tabs, qvap, qcld )

  tabs = tabs + d_ta
  qvap = qvap + d_qv
  qcld = qcld + d_qc

!-----------------------------------------------------------------------
! get changes in ta, qv, qc and qr due to autoconversion
!-----------------------------------------------------------------------

  call get_qsat ( pres, tabs, qsat, wice )

  call evaporate ( Mgrid, pres, tabs, qvap, qcld, dpm, rain, delt )

  tabs = tabs + d_ta
  qvap = qvap + d_qv
  qcld = qcld + d_qc

!-----------------------------------------------------------------------
! precipitation fields
!-----------------------------------------------------------------------

  if ( do_ice ) then
     do j=jbc,jec
        do i=ibc,iec
           wsleet = tanh((tsleet - tabs(i,j,ked))/(tsleet - tfreeze))
           snow(i,j) = rain(i,j)*(0.5 + 0.5*wsleet)
           rain(i,j) = rain(i,j)*(0.5 - 0.5*wsleet)
        enddo
     enddo
  endif

!-----------------------------------------------------------------------
! write diagnostic fields
!-----------------------------------------------------------------------

  d_ta = (tabs - tabs1)/delt
  d_qv = (qvap - qvap1)/delt
  d_qc = (qcld - qcld1)/delt
  rain = rain/delt
  snow = snow/delt

  call ncdf_fms_write (Mgrid, Time, id_rain_mp, gridm, rain)
  call ncdf_fms_write (Mgrid, Time, id_snow_mp, gridm, snow)
!!$  call ncdf_fms_write (Mgrid, Time, id_tdt_mp, gridm, d_ta)
!!$  call ncdf_fms_write (Mgrid, Time, id_qdt_mp, gridm, d_qv)
!!$  call ncdf_fms_write (Mgrid, Time, id_ldt_mp, gridm, d_qc)

  return
end subroutine kessler_mp

!#######################################################################

subroutine condense ( Mgrid, pp, ta, qv, qc )

type (horiz_grid_type),                           intent(in)     ::    &
                                                                Mgrid

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),             intent (in)    ::    &
                                                       pp, ta, qv, qc
        
!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

real :: cond, qdif, heat
real, parameter :: qvfrac=0.2
integer :: i, j, k
  
!-----------------------------------------------------------------------
! compute condensation [dqdT = d(qsat)/dT]
!-----------------------------------------------------------------------

  call get_dlnesatdt ( ta, dqdt, wice )
  dqdt = qsat * dqdt * (1. + qsat/eps)

  do k=kbd+1,ked
     do j=jbc,jec
        do i=ibc,iec
           qdif = qv(i,j,k) - qsat(i,j,k)            ! supersaturation
           heat = hlat(i,j,k) / cp_air
           cond = max(qdif/(1.0 + dqdt(i,j,k)*heat),0.)
!!           cond = max(qdif/(1.0 + dqdt(i,j,k)*heat),                   &
!!              -min(0.,max(qc(i,j,k),-qv(i,j,k)*qvfrac))) ! fill up qc<0
           d_qv(i,j,k) = -cond
           d_qc(i,j,k) =  cond
           d_ta(i,j,k) =  cond*heat
        enddo
     enddo
  enddo

  return
end subroutine condense

!#######################################################################

subroutine evaporate ( Mgrid, pp, ta, qv, qc, dp, rain, dt )

!-----------------------------------------------------------------------
! rainwater parameterizations
!-----------------------------------------------------------------------

type (horiz_grid_type),                           intent(in)     ::    &
                                                                Mgrid

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),             intent (in)    ::    &
                                                   pp, ta, qv, qc, dp
        
real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed),             &
                                                  intent (out)   ::    &
                                                                 rain
        
real,                                             intent (in)    ::    &
                                                                   dt
    
!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------
  
real, parameter :: q0=0.01
real, parameter :: cthrice=0.0e-03, ricefac=0.05

real :: rliq, rice, qdif, heat
real :: accr, evap, evapdt, autodt
real :: qv1, qc1, qr1

integer :: i, j, k

!-----------------------------------------------------------------------
! compute autoconversion of cloud water
!-----------------------------------------------------------------------

  rliq = rauto
  rice = rauto*ricefac

  do k=kbd+1,ked
     do j=jbc,jec
        do i=ibc,iec
           qc1 = qc(i,j,k)
           auto(i,j,k) = rliq*dim(qc1*(1.0 - wice(i,j,k)), cthresh)    &
                       + rice*dim(qc1*wice(i,j,k)        , cthrice)
        enddo
     enddo
  enddo

  if ( do_evap ) then

!-----------------------------------------------------------------------
! include accretion and re-evaporation of rain ...
!-----------------------------------------------------------------------

     call get_dlnesatdt ( ta, dqdt, wice )
     dqdt = qsat * dqdt * (1. + qsat/eps)

     rain = 0.
     do k=kbd+1,ked
        do j=jbc,jec
           do i=ibc,iec
              qv1 = qv(i,j,k)
              qc1 = qc(i,j,k)
              qr1  = dim(rain(i,j)/dp(i,j,k), rthresh)/q0 + tiny
              accr  = dim(qc1, 0.0)*raccr * qr1**eaccr
              qdif  = qsat(i,j,k) - qv1                ! subsaturation

              evap  = revap * qdif*(q0/qsat(i,j,k)) * qr1**eevap
              heat  = hlat(i,j,k) / cp_air
              qdif = qdif/(1.0 + dqdt(i,j,k)*heat)

              autodt = max(0.0, min(0.9*qc1, (accr + auto(i,j,k))*dt))
              evapdt = max(0.0, min(0.9*qdif, 0.9*qr1, evap*dt))

              d_qv(i,j,k) = evapdt
              d_qc(i,j,k) = -autodt
              d_ta(i,j,k) = -evapdt*heat
              rain(i,j) = rain(i,j) + (autodt - evapdt)*dp(i,j,k)
           enddo
        enddo
     enddo

  else if ( do_auto ) then

!-----------------------------------------------------------------------
! ... or use autoconversion only
!-----------------------------------------------------------------------

     rain = 0.
     do k=kbd+1,ked
        do j=jbc,jec
           do i=ibc,iec
              qc1 = qc(i,j,k)
              autodt = max(0.0, min(0.9*qc1, auto(i,j,k)*dt))
              d_qc(i,j,k) = -autodt
              rain(i,j) = rain(i,j) + autodt*dp(i,j,k)
           enddo
        enddo
     enddo
     d_qv = 0.
     d_ta = 0.

  else

     rain = 0.
     do k=kbd+1,ked
        do j=jbc,jec
           do i=ibc,iec
              d_qc(i,j,k) = -qc(i,j,k)
              rain(i,j) = rain(i,j) + qc(i,j,k)*dp(i,j,k)
           enddo
        enddo
     enddo
     d_qv = 0.
     d_ta = 0.

  endif

  rain = rain/Grav

  return
end subroutine evaporate

!#######################################################################

subroutine kessler_init ( Mgrid, Time )

type (horiz_grid_type),                     intent (in)    ::          &
                                                                Mgrid
       
type (time_type),                           intent (in)    ::          &
                                                                 Time

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

integer :: io, ierr
integer, dimension(2) :: axid_mm
integer, dimension(3) :: axid_mmm

character (len=16) :: name, units
character (len=64) :: longname

namelist /kessler_nml/ revap, raccr, rauto, cthresh, rthresh,          &
                       do_evap, do_auto, do_ice

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

integer :: ibd, ied
integer :: jbd, jed
integer :: n, nt

  ibd = Mgrid%ibd
  ied = Mgrid%ied
  jbd = Mgrid%jbd
  jed = Mgrid%jed
  kbd = Mgrid%kbd
  ked = Mgrid%ked

  ibc = Mgrid%ibc
  iec = Mgrid%iec
  jbc = Mgrid%jbc
  jec = Mgrid%jec

!-----------------------------------------------------------------------
! read namelist
!-----------------------------------------------------------------------

  read (input_nml_file, nml=kessler_nml, iostat=io)
  ierr = check_nml_error(io,'kessler_nml')
 
!-----------------------------------------------------------------------
! write version number and namelist to logfile
!-----------------------------------------------------------------------

  call write_version_number (version, tag)
  if (mpp_pe() == mpp_root_pe())                                       &
                                      write (stdlog(), nml=kessler_nml)

  allocate ( qsat(ibd:ied, jbd:jed, kbd:ked) )
  allocate ( auto(ibd:ied, jbd:jed, kbd:ked) )
  allocate ( wice(ibd:ied, jbd:jed, kbd:ked) )
  allocate ( hlat(ibd:ied, jbd:jed, kbd:ked) )
  allocate ( dqdt(ibd:ied, jbd:jed, kbd:ked) )

  allocate ( d_ta(ibd:ied, jbd:jed, kbd:ked) )
  allocate ( d_qv(ibd:ied, jbd:jed, kbd:ked) )
  allocate ( d_qc(ibd:ied, jbd:jed, kbd:ked) )

  allocate ( tabs1(ibd:ied, jbd:jed, kbd:ked) )
  allocate ( qvap1(ibd:ied, jbd:jed, kbd:ked) )
  allocate ( qcld1(ibd:ied, jbd:jed, kbd:ked) )

!-----------------------------------------------------------------------
! register diagnostics
!-----------------------------------------------------------------------

  axid_mmm = (/axid(xm_),axid(ym_),axid(zm_)/)
  axid_mm  = (/axid(xm_),axid(ym_)/)

  name = "rain_mp"
  longname = "rain rate from Kessler"
  units = 'kg/m2/s'
  id_rain_mp = register_diag_field (                                   &
                      model, name, axid_mm, Time,                      &
                      longname, units, silly, range=(/-10., 10./) )

  name = "snow_mp"
  longname = "snow rate from Kessler"
  units = 'kg/m2/s'
  id_snow_mp = register_diag_field (                                   &
                      model, name, axid_mm, Time,                      &
                      longname, units, silly, range=(/-10., 10./) )

!!$  name = "tdt_mp"
!!$  longname = "temperature tendency from Kessler"
!!$  units = 'K/s'
!!$  id_tdt_mp = register_diag_field (                                    &
!!$                      model, name, axid_mmm, Time,                     &
!!$                      longname, units, silly, range=(/-1.e4, 1.e4/) )
!!$
!!$  name = "qdt_mp"
!!$  longname = "water vapor tendency from Kessler"
!!$  units = 'kg/kg/s'
!!$  id_qdt_mp = register_diag_field (                                    &
!!$                      model, name, axid_mmm, Time,                     &
!!$                      longname, units, silly, range=(/-1.e-3, 1.e-3/) )
!!$
!!$  name = "ldt_mp"
!!$  longname = "condensate tendency from Kessler"
!!$  units = 'kg/kg/s'
!!$  id_ldt_mp = register_diag_field (                                    &
!!$                      model, name, axid_mmm, Time,                     &
!!$                      longname, units, silly, range=(/-1.e-5, 1.e-5/) )

  return
end subroutine kessler_init

!#######################################################################

subroutine kessler_end ( Mgrid, time, ltime )

type (horiz_grid_type)                                       ::        &
                                                                Mgrid

real,                                         intent (in)    ::        &
                                                                 time

integer,                                      intent (in)    ::        &
                                                                ltime

  return
end subroutine kessler_end

!#######################################################################
   
subroutine error_handler ( message ) 
character(len=*), intent(in) :: message
   
   call error_mesg (module, message, FATAL)

end subroutine error_handler

!#######################################################################

end module zetac_kessler_mod
