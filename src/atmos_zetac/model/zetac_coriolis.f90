module zetac_coriolis_mod

!-----------------------------------------------------------------------
! module information
!-----------------------------------------------------------------------

use mpp_mod,                     only : mpp_pe, mpp_root_pe
use fms_mod,                     only : write_version_number

use zetac_coriolis_param_mod,    only : coriolis_param
use zetac_horiz_grid_type_mod,   only : horiz_grid_type
use zetac_horiz_metric_type_mod, only : horiz_metric_type
use zetac_time_pointers_mod,     only : ntime

!-----------------------------------------------------------------------
! public interfaces
!-----------------------------------------------------------------------

implicit none
private
public coriolis_at_u, coriolis_at_v, coriolis_init
                                     
character(len=*), parameter :: module='zetac_coriolis_mod'

real, save, allocatable, dimension(:) :: ffu, cormetu
real, save, allocatable, dimension(:) :: ffv, cormetv
real, save, allocatable, dimension(:,:,:) :: vden

integer :: ibc, iec
integer :: jbc, jec
integer :: kbc, kec

!-----------------------------------------------------------------------
! version control information
!-----------------------------------------------------------------------

character(len=128)  :: version =  '$Id: zetac_coriolis.f90,v 1.1.2.2.2.3 2005/06/16 18:07:28 stg Exp $'
character(len=128)  :: tag     =  '$Name:  $'

contains

!#######################################################################

subroutine coriolis_at_u ( Mgrid, denu, denv, uu, vv, coriolis )

!-----------------------------------------------------------------------
! calculate zonal velocity tendency due to coriolis force
!-----------------------------------------------------------------------
 
type (horiz_grid_type),                     intent (in)    ::          &
                                                                Mgrid

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),       intent (in)    ::          &
                                                   denu, denv, uu, vv

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),       intent (out)   ::          &
                                                             coriolis

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

real    :: vatu
integer :: i, j, k

  do k=kbc,kec
     do j=jbc-1,jec
        do i=ibc,iec+1
           vden(i,j,k) = vv(i,j,k)*sqrt(denv(i,j,k))
        enddo
     enddo
  enddo
  do k=kbc,kec
     do j=jbc,jec
        do i=ibc,iec
           vatu = 0.25*(vden(i,j  ,k) + vden(i+1,j  ,k)                &
                      + vden(i,j-1,k) + vden(i+1,j-1,k))/              &
                                   sqrt(denu(i,j,k))
           coriolis(i,j,k) = (ffu(j) + cormetu(j)*uu(i,j,k))*vatu
        enddo
     enddo
  enddo

  return
end subroutine coriolis_at_u

!#######################################################################

subroutine coriolis_at_v ( Mgrid, denu, denv, uu, vv, coriolis )

!-----------------------------------------------------------------------
! calculate meridional velocity tendency due to coriolis force
!-----------------------------------------------------------------------
  
type (horiz_grid_type),                     intent (in)    ::          &
                                                                Mgrid

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),       intent (in)    ::          &
                                                   denu, denv, uu, vv 

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),       intent (out)   ::          &
                                                             coriolis

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

real    :: uatv
integer :: i, j, k

  do k=kbc,kec
     do j=jbc,jec+1
        do i=ibc-1,iec
           vden(i,j,k) = uu(i,j,k)*sqrt(denu(i,j,k))
        enddo
     enddo
  enddo
  do k=kbc,kec
     do j=jbc,jec
        do i=ibc,iec
           uatv = .25*(vden(i,j  ,k) + vden(i-1,j  ,k)                 &
                     + vden(i,j+1,k) + vden(i-1,j+1,k))/               &
                                  sqrt(denv(i,j,k))
           coriolis(i,j,k) = -(ffv(j) + cormetv(j)*uatv)*uatv
        enddo
     enddo
  enddo
  
  return
end subroutine coriolis_at_v

!#######################################################################

subroutine coriolis_init ( Mgrid, Hmetric )

type (horiz_grid_type),                     intent (in)    ::          &
                                                                Mgrid

type (horiz_metric_type),                   intent (in)    ::          &
                                                              Hmetric

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

real    :: rlatmin, rlatmax

real, dimension(Mgrid%jbg:Mgrid%jeg) :: ulat
real, dimension(Mgrid%jbg:Mgrid%jeg) :: vlat

!-----------------------------------------------------------------------
! write version number to logfile
!-----------------------------------------------------------------------

  call write_version_number (version, tag)

  ibc = Mgrid%ibc
  iec = Mgrid%iec
  jbc = Mgrid%jbc
  jec = Mgrid%jec
  kbc = Mgrid%kbc
  kec = Mgrid%kec

  rlatmin = Hmetric%rlatmin
  rlatmax = Hmetric%rlatmax
  ulat = Hmetric%ulat
  vlat = Hmetric%vlat

  allocate (    ffu(Mgrid%jbg:Mgrid%jeg))
  allocate (    ffv(Mgrid%jbg:Mgrid%jeg))
  allocate (cormetu(Mgrid%jbg:Mgrid%jeg))
  allocate (cormetv(Mgrid%jbg:Mgrid%jeg))

  allocate (vden(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed, Mgrid%iz))

  call coriolis_param ( Mgrid, Hmetric, ulat, ffu, cormetu )
  call coriolis_param ( Mgrid, Hmetric, vlat, ffv, cormetv )

  return
end subroutine coriolis_init

!#######################################################################

end module zetac_coriolis_mod

