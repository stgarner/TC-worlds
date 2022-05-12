module zetac_write_coldstart_mod

!-----------------------------------------------------------------------
! module information
!-----------------------------------------------------------------------

use mpp_mod,    only : mpp_pe, mpp_root_pe
use fms_mod,    only : file_exist, close_file, write_version_number
use mpp_io_mod, only : fieldtype, axistype, mpp_get_id,                &
                       mpp_write, mpp_write_meta, mpp_open, mpp_close, &
                       MPP_NETCDF, MPP_SINGLE, MPP_OVERWR
use mpp_domains_mod,             only : domain1D,                      &
                                        mpp_get_domain_components

use constants_mod,               only : kappa
use diag_manager_mod,            only : get_base_date
use field_manager_mod,           only : model_atmos
use time_manager_mod,            only : time_type

use zetac_axis_names_mod,        only : xm_, ym_
use zetac_horiz_grid_type_mod,   only : horiz_grid_type
use zetac_horiz_metric_type_mod, only : horiz_metric_type
use zetac_ncdf_io_mod,           only : ncdf_io_init
use zetac_restart_mod,           only : write_restart,                 &
                                        write_restart_init
use zetac_spec_hum_mod,          only : get_spec_hum
use zetac_update_halos_mod,      only : update_halos
use zetac_vert_metric_type_mod,  only : vert_metric_type

!-----------------------------------------------------------------------
! public interfaces
!-----------------------------------------------------------------------

implicit none
private
public write_coldstart, coldstart_init, coldstart_end

character(len=*), parameter :: module='zetac_write_coldstart_mod'
character(len=*), parameter :: filename_sst = 'OUTPUT/sst'

logical :: found_landmask

real, allocatable, dimension(:,:)   :: ppref
real, allocatable, dimension(:,:)   :: qvref
real, allocatable, dimension(:,:,:) :: qvpert
real, allocatable, dimension(:,:,:) :: qcpert
real, allocatable, dimension(:,:,:) :: pppert
real, allocatable, dimension(:,:,:) :: oopert
real, allocatable, dimension(:,:,:) :: dppert

real, allocatable, dimension(:) :: zetam, zetaw

type(time_type), save :: Time_step
type(domain1D),  save :: Domain_x, Domain_y
type(axistype),  save :: Axis_x, Axis_y, Axis_t
type(fieldtype), save :: Field_sst

integer :: unit_mask, unit_sst

integer :: ibc, iec, ibh, ieh
integer :: jbc, jec, jbh, jeh
integer :: kbd, ked

real :: ptop

!-----------------------------------------------------------------------
! version control information
!-----------------------------------------------------------------------

character(len=128)  :: version =  '$Id: zetac_write_coldstart.f90,v 1.1.2.4.2.8 2005/07/05 02:56:53 stg Exp $'
character(len=128)  :: tag     =  '$Name:  $'

contains

!#######################################################################

subroutine write_coldstart (                                           &
                                                 Mgrid, Vmetric, Time, &
                                                           topog, sst, &
                                    psref, tsref, uuref, taref, rhref, &
                               pspert, uupert, vvpert, tapert, rhpert )

!-----------------------------------------------------------------------
! writes out fixed and initial time-dependent arrays to restart file
!-----------------------------------------------------------------------

type (horiz_grid_type)                                       ::        &
                                                                Mgrid

type (vert_metric_type),                      intent (in)    ::        &
                                                              Vmetric

type (time_type),                             intent (in)    ::        &
                                                                 Time

real,dimension(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed),               &
                                              intent (inout) ::        &
                                                           topog, sst

real,dimension(Mgrid%jbg:Mgrid%jeg),          intent (inout) ::        &
                                                         psref, tsref

real,dimension(Mgrid%jbg:Mgrid%jeg,Mgrid%iz), intent (inout) ::        &
                                                  uuref, taref, rhref

real,dimension(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed),               &
                                              intent (inout) ::        &
                                                               pspert

real,dimension(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed,                &
               Mgrid%kbd:Mgrid%ked),          intent (inout) ::        &
                                       uupert, vvpert, tapert, rhpert

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

integer :: i, j, k

!-----------------------------------------------------------------------
! combine perturbation with reference atmosphere 
!-----------------------------------------------------------------------

  call combine_fields (                                                &
                                                       Mgrid, Vmetric, &
                                           psref, uuref, taref, rhref, &
                                       pspert, uupert, tapert, rhpert )

!-----------------------------------------------------------------------
! get specific humidity from relative humidity
!-----------------------------------------------------------------------

  call get_spec_hum (                                                  &
                                                       Mgrid, Vmetric, &
                                           ppref, taref, rhref, qvref )

  call get_spec_hum (                                                  &
                                                       Mgrid, Vmetric, &
                                       pppert, tapert, rhpert, qvpert )

!-----------------------------------------------------------------------
! create cold-start file
!-----------------------------------------------------------------------

  call update_halos ( Mgrid, pspert )
  call update_halos ( Mgrid, pppert )
  call update_halos ( Mgrid, uupert )
  call update_halos ( Mgrid, vvpert )
  call update_halos ( Mgrid, tapert )
  call update_halos ( Mgrid, qvpert )

  do k=kbd+1,ked                   ! dry mass in layers
     dppert(:,:,k) = (pppert(:,:,k) - pppert(:,:,k-1))/                &
                                (1. + qvpert(:,:,k))
  enddo
  dppert(:,:,kbd) = dppert(:,:,kbd+1)

  do j=Mgrid%jbd,Mgrid%jed         ! overwrite sst
     sst(:,j) = tsref(j)
  enddo

  oopert = 0.
  qcpert = 0.

  call write_restart (                                                 &
                                                    Mgrid, Time, Time, &
                                                           topog, sst, &
                                           ppref, uuref, taref, qvref, &
       pspert, dppert, uupert, vvpert, oopert, tapert, qvpert, qcpert )

  return
end subroutine write_coldstart

!#######################################################################

subroutine combine_fields (                                            &
                                                       Mgrid, Vmetric, &
                                           psref, uuref, taref, rhref, &
                                       pspert, uupert, tapert, rhpert )

!-----------------------------------------------------------------------
! combines reference state with perturbation
!-----------------------------------------------------------------------

type (horiz_grid_type),                        intent (in)      ::     &
                                                                Mgrid

type (vert_metric_type),                      intent (in)    ::        &
                                                              Vmetric

real, dimension(Mgrid%jbg:Mgrid%jeg),                                  &
                                               intent (in)      ::     &
                                                                psref

real, dimension(Mgrid%jbg:Mgrid%jeg, Mgrid%iz),                        &
                                               intent (in)      ::     &
                                                  uuref, taref, rhref

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed),             &
                                               intent (inout)   ::     &
                                                               pspert

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),          intent (inout)   ::     &
                                               uupert, tapert, rhpert

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

integer :: i, j, k

  do k=1,Mgrid%iz
     ppref(:,k) = ptop + (psref(:) - ptop)*Vmetric%zetaw(k)
  enddo

  do j=jbh,jeh
     do i=ibh,ieh
        pspert(i,j) = pspert(i,j) + psref(j)
     enddo
  enddo

  do k=1,Mgrid%iz
     pppert(:,:,k) = ptop + (pspert(:,:) - ptop)*Vmetric%zetaw(k)
  enddo

  do k=kbd,ked
     do j=jbh,jeh
        do i=ibh,ieh
           uupert(i,j,k) = uupert(i,j,k) + uuref(j,k)
           tapert(i,j,k) = tapert(i,j,k) + taref(j,k)
           rhpert(i,j,k) = rhpert(i,j,k) + rhref(j,k)
        enddo
     enddo
  enddo

  return
end subroutine combine_fields

!#######################################################################

subroutine coldstart_init (                                            &
                                              Mgrid, Hmetric, Vmetric, &
                                   rlonmin, rlonmax, rlatmin, rlatmax, &
                                         do_cartesian, do_cylindrical )

type (horiz_grid_type),                            intent (in)    ::   &
                                                                Mgrid

type (horiz_metric_type),                          intent (in)    ::   &
                                                              Hmetric

type (vert_metric_type),                           intent (in)    ::   &
                                                              Vmetric

logical,                                           intent (in)    ::   &
                                         do_cartesian, do_cylindrical

real,                                              intent (in)    ::   &
                                   rlonmin, rlonmax, rlatmin, rlatmax

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

character(len=2)  :: mstring, dstring
character(len=34) :: string
character(len=64) :: calendar_type, name, units
integer :: length, id
integer :: year, month, day, date(6)

!-----------------------------------------------------------------------
! write version number to logfile
!-----------------------------------------------------------------------

  call write_version_number (version, tag)

!-----------------------------------------------------------------------
! initialize ncdf_io
!-----------------------------------------------------------------------

  call ncdf_io_init ( Mgrid, Time_step )

!-----------------------------------------------------------------------
! initialize restart routine
!-----------------------------------------------------------------------

  call write_restart_init ( Mgrid, Hmetric, Vmetric )

!-----------------------------------------------------------------------
! allocate work arrays
!-----------------------------------------------------------------------

  allocate (ppref(Mgrid%jbg:Mgrid%jeg,Mgrid%iz))
  allocate (qvref(Mgrid%jbg:Mgrid%jeg,Mgrid%iz))

  allocate (qvpert(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed,Mgrid%iz))
  allocate (pppert(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed,Mgrid%iz))
  allocate (qcpert(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed,Mgrid%iz))
  allocate (oopert(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed,Mgrid%iz))
  allocate (dppert(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed,Mgrid%iz))

  ibc = Mgrid%ibc
  iec = Mgrid%iec
  ibh = Mgrid%ibh
  ieh = Mgrid%ieh

  jbc = Mgrid%jbc
  jec = Mgrid%jec
  jbh = Mgrid%jbh
  jeh = Mgrid%jeh

  kbd = Mgrid%kbd
  ked = Mgrid%ked

  ptop = Vmetric%ptop

  allocate (zetam(Mgrid%iz))
  allocate (zetaw(Mgrid%iz))

  zetam = Vmetric%zetam
  zetaw = Vmetric%zetaw

  return
end subroutine coldstart_init

!#######################################################################

subroutine coldstart_end

  return
end subroutine coldstart_end

!#######################################################################

end module zetac_write_coldstart_mod
