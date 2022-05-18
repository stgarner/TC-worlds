module zetac_ls_forcing_mod

!-----------------------------------------------------------------------
!
!       This module controls the following functions:
!
!       (1) Simple cooling (Combination of Newtonian cooling
!                    and prescribed cooling rates
!
!       (2) idealized large-scale forcing 
!
!-----------------------------------------------------------------------

use zetac_horiz_grid_type_mod, only : horiz_grid_type
use fms_mod,          only:  mpp_pe, mpp_root_pe, file_exist,          &
                             error_mesg, check_nml_error, FATAL,       &
                             open_file, close_file,                    &
                             open_namelist_file, stdlog,               &
                             write_version_number
use mpp_mod,          only:  input_nml_file
use constants_mod,    only:  cp_air, grav, hlv
use time_manager_mod, only : time_type
use diag_manager_mod, only : register_diag_field

use zetac_axes_mod,       only : axid, gridm
use zetac_axis_names_mod, only : xm_, ym_, zm_
use zetac_ncdf_io_mod,    only : ncdf_fms_write
use zetac_vert_metric_type_mod,  only : vert_metric_type

implicit none
private

public ls_forcing, ls_forcing_init

character (len=*),  parameter :: module='zetac_ls_forcing_mod'
character (len=16), parameter :: model='atmos_mod'
character (len=16) :: name, units
character (len=64) :: longname

real, parameter :: silly=-9999.9

! namelist:
real :: strat_temp = 210.
real ::  crit_temp = 250.
real :: nudge_time = 1.0
real ::  cool_time = 1.5

namelist /ls_forcing_nml/ strat_temp, crit_temp, cool_time, nudge_time

integer ::  k, iz
real :: missing_value = -999.

integer :: id_tdt_ls

character(len=16) :: mod_name = 'zetac_ls_forcing_mod'

real, allocatable, dimension(:) :: nudge_tend, cool_tend, temp_strat
real, allocatable, dimension(:) :: zetam

character(len=128) :: version = '$Id: zetac_ls_forcing.f90,v 1.1.2.1.2.2.2.2 2004/09/15 16:03:45 ck Exp $'
character(len=128) :: tagname = '$Name:  $'

contains

!#######################################################################

subroutine ls_forcing (                                                &
                                               Mgrid, Time, Time_next, &
                                                               ta, qv, &
                                                                 delt )

type (horiz_grid_type),                     intent (in)    ::          &
                                                                Mgrid

type (time_type),                           intent (in)    ::          &
                                                      Time, Time_next

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),       intent(inout)  ::          &
                                                               ta, qv

real,                                       intent(in)     ::    delt

!-----------------------------------------------------------------------
! local allocations  
!-----------------------------------------------------------------------

real, dimension(size(ta,1),size(ta,2),size(ta,3)) :: tdt, qdt

  call simple_cooling ( delt, ta, tdt )
  ta = ta + tdt*delt

!-----------------------------------------------------------------------
! write diagnostic fields
!-----------------------------------------------------------------------

  call ncdf_fms_write (Mgrid, Time, id_tdt_ls, gridm, tdt)

  return
 end subroutine ls_forcing

!#######################################################################

 subroutine simple_cooling ( dt, temp, dtdt )

  real, intent(in) :: dt
  real, intent(in),  dimension(:,:,:)  :: temp
  real, intent(out), dimension(:,:,:)  :: dtdt

!-----------------------------------------------------------------------
! simple cooling
!-----------------------------------------------------------------------

!!$  do k=1,iz
!!$     WHERE ( temp(:,:,k) > 285. )
!!$        dtdt(:,:,k) = cool_tend(k)!*0.5
!!$     ELSEWHERE ( temp(:,:,k) > crit_temp )
!!$        dtdt(:,:,k) = cool_tend(k)
!!$     ELSEWHERE
!!$        dtdt(:,:,k) = nudge_tend(k)*(temp_strat(k) - temp(:,:,k))/     &
!!$                (1. + nudge_tend(k)*dt)
!!$     ENDWHERE
!!$  end do
  do k=1,iz
     if ( zetam(k) > 0.25 ) then
        dtdt(:,:,k) = cool_tend(k)
     else
        dtdt(:,:,k) = nudge_tend(k)*(temp_strat(k) - temp(:,:,k))/     &
                (1. + nudge_tend(k)*dt)
     endif
  end do

 end subroutine simple_cooling

!#######################################################################

 subroutine ls_forcing_init ( Mgrid, Vmetric, Time )

type (horiz_grid_type),                     intent (in)    ::          &
                                                                Mgrid
       
type (vert_metric_type),                    intent (in)    ::          &
                                                              Vmetric

type (time_type),                           intent (in)    ::          &
                                                                 Time

integer :: ierr, io
integer, dimension(2) :: axid_mm
integer, dimension(3) :: axid_mmm

!-----------------------------------------------------------------------
! namelist (read & write) 
!-----------------------------------------------------------------------

  read (input_nml_file, nml=ls_forcing_nml, iostat=io)
  ierr = check_nml_error(io,'ls_forcing_nml')

!-----------------------------------------------------------------------
! write version number and namelist to logfile
!-----------------------------------------------------------------------

  call write_version_number (version, tagname)
  if (mpp_pe() == mpp_root_pe())                                       &
                                   write (stdlog(), nml=ls_forcing_nml)

  iz = Mgrid%iz

  !-----------------------------------------------------------------------
!  simple cooling
!-----------------------------------------------------------------------

  allocate(nudge_tend(iz))
  allocate(cool_tend(iz))
  allocate(temp_strat(iz)) 

  do k=1,iz
     nudge_tend(k) = 1./(nudge_time*86400.)
     cool_tend(k) = -1./(cool_time*86400.)
     temp_strat(k) = strat_temp
  enddo

!-----------------------------------------------------------------------
! register diagnostics
!-----------------------------------------------------------------------

  axid_mmm = (/axid(xm_),axid(ym_),axid(zm_)/)
  axid_mm  = (/axid(xm_),axid(ym_)/)

  name = "tdt_ls"
  longname = "cooling rate from ls_forcing"
  units = 'K/s'
  id_tdt_ls = register_diag_field (                                    &
                      model, name, axid_mmm, Time,                     &
                      longname, units, silly, range=(/-1., 1./) )

  allocate (zetam(iz))
  zetam = Vmetric%zetam

  return
end subroutine ls_forcing_init

!#######################################################################

end module zetac_ls_forcing_mod

