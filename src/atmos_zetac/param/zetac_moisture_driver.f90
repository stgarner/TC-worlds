module zetac_moisture_driver_mod

!-----------------------------------------------------------------------
! module information
!-----------------------------------------------------------------------

use mpp_mod,    only : mpp_pe, mpp_root_pe, mpp_broadcast,             &
                       mpp_clock_id, mpp_clock_begin,                  &
                       mpp_clock_end, input_nml_file
use fms_mod,    only : file_exist, write_version_number,               &
                       close_file, open_namelist_file,                 &
                       check_nml_error, stdlog, stdout

use diag_manager_mod,            only : register_diag_field
use field_manager_mod,           only : model_atmos
use time_manager_mod,            only : time_type, get_time,           &
                                        operator(-)
use constants_mod,               only : grav

use zetac_calc_var_tend_mod,     only : get_gz
use zetac_convert_var_mod,       only : get_tv, get_pp
use zetac_extrap_var_mod,        only : extrap_tracer
use zetac_horiz_grid_type_mod,   only : horiz_grid_type
use zetac_horiz_metric_type_mod, only : horiz_metric_type
use zetac_kessler_mod,           only : kessler_mp,                    &
                                        kessler_init, kessler_end
use zetac_ras_mod,               only : ras, ras_init, ras_end
use zetac_ncdf_io_mod,           only : ncdf_fms_write,                &
                                        ncdf_read, ncdf_write
use zetac_update_halos_mod,      only : update_halos
use zetac_vert_metric_type_mod,  only : vert_metric_type

!-----------------------------------------------------------------------
! public interfaces
!-----------------------------------------------------------------------

implicit none
private
public moisture_driver,                                                &
       moisture_driver_init, moisture_driver_end

integer :: zetac_moisture_driver_clock

type(time_type), save :: Time_base

character (len=*),  parameter :: module='zetac_moisture_driver_mod'
character (len=16), parameter :: model='atmos_mod'

real, parameter :: silly=-99999.9

integer :: kbd, ked

!-----------------------------------------------------------------------
! namelist parameters
!-----------------------------------------------------------------------

logical :: do_kessler=.true.
logical :: do_ras=.false.
logical :: do_CC
real :: qvmin=0.

!-----------------------------------------------------------------------
! version control information
!-----------------------------------------------------------------------

character(len=128)  :: version =  '$Id: zetac_moisture_driver.f90,v 1.1.2.10.2.6.2.7 2005/07/06 19:53:13 stg Exp $'
character(len=128)  :: tag     =  '$Name:  $'

contains

!#######################################################################

subroutine moisture_driver (                                           &
                                              Mgrid, Hmetric, Vmetric, &
                                                                 Time, &
                                               psfc, tabs, qvap, qcld, &
                                                      dpm, rain, snow, &
                                                                 delt )

type (horiz_grid_type),                       intent (in)    ::        &
                                                                Mgrid

type (horiz_metric_type),                     intent (in)    ::        &
                                                              Hmetric
       
type (vert_metric_type),                      intent (in)    ::        &
                                                              Vmetric
       
type (time_type),                             intent (in)    ::        &
                                                                 Time

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed),             &
                                              intent (in)    ::        &
                                                                 psfc

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),         intent (inout) ::        &
                                                     tabs, qvap, qcld

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),         intent (in)    ::        &
                                                                  dpm

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed),             &
                                              intent (out)   ::        &
                                                           rain, snow

real,                                         intent (in)    ::        &
                                                                 delt

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked)                  ::  pfull, phalf
real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked)                  ::  zfull, zhalf

!-----------------------------------------------------------------------
! density field
!-----------------------------------------------------------------------

  call mpp_clock_begin ( zetac_moisture_driver_clock )

  call get_pp ( psfc, Vmetric%zetam, Vmetric%ptop, pfull )
  call get_pp ( psfc, Vmetric%zetaw, Vmetric%ptop, phalf )

  call get_gz ( Mgrid, phalf, pfull, tabs, qvap, zhalf, zfull )

  zhalf = zhalf/Grav
  zfull = zfull/Grav

!-----------------------------------------------------------------------
! "relaxed Arakawa-Shubert" convection
!-----------------------------------------------------------------------

  if (do_ras) then
     call ras (                                           Mgrid, Time, &
                                           zfull, pfull, zhalf, phalf, &
                                                     tabs, qvap, qcld, &
                                                                 delt )
  endif

!-----------------------------------------------------------------------
! kessler microphysics
!-----------------------------------------------------------------------

  if (do_kessler) then
     call kessler_mp (                                                 &
                                                          Mgrid, Time, &
                                                           pfull, dpm, &
                                                     tabs, qvap, qcld, &
                                                           rain, snow, &
                                                                 delt )
  endif
  
  call mpp_clock_end ( zetac_moisture_driver_clock )

  return
end subroutine moisture_driver

!#######################################################################

subroutine moisture_driver_init ( Mgrid, Time, qvmin_out )

type (horiz_grid_type),                      intent (in)    ::         &
                                                                Mgrid

type (time_type),                            intent (in)    ::         &
                                                                 Time

real,                                        intent (out)   ::         &
                                                            qvmin_out

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

integer :: io, ierr

namelist /moisture_nml/ do_kessler, do_ras, qvmin, do_CC

  kbd = Mgrid%kbd
  ked = Mgrid%ked

!-----------------------------------------------------------------------
! read namelist
!-----------------------------------------------------------------------

  read(input_nml_file, nml=moisture_nml, iostat=io)
  ierr = check_nml_error(io, 'moisture_nml')

  zetac_moisture_driver_clock = mpp_clock_id( 'zetac_moisture' )

!-----------------------------------------------------------------------
! write version number and namelist to logfile
!-----------------------------------------------------------------------

  call write_version_number (version, tag)
  if (mpp_pe() == mpp_root_pe()) write (stdlog(), nml=moisture_nml)

  qvmin_out = qvmin

!-----------------------------------------------------------------------
! initialize RAS
!-----------------------------------------------------------------------

  if (do_ras) call ras_init ( Mgrid, Time )

!-----------------------------------------------------------------------
! initialize Kessler microphysics
!-----------------------------------------------------------------------

  if (do_kessler) call kessler_init ( Mgrid, Time )

  return
end subroutine moisture_driver_init

!#######################################################################

subroutine moisture_driver_end ( Mgrid, time, ltime )

type (horiz_grid_type)                                      ::         &
                                                                Mgrid

real,                                       intent (in)     ::         &
                                                                 time

integer,                                    intent (in)     ::         &
                                                                ltime

!-----------------------------------------------------------------------
! finalize RAS
!-----------------------------------------------------------------------

  if (do_ras) call ras_end ( Mgrid, time, ltime )

!-----------------------------------------------------------------------
! finalize microphysics
!-----------------------------------------------------------------------

  if (do_kessler) call kessler_end ( Mgrid, time, ltime )

  return
end subroutine moisture_driver_end

!#######################################################################

end module zetac_moisture_driver_mod
