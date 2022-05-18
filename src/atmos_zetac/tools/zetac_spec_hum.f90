module zetac_spec_hum_mod

!-----------------------------------------------------------------------
! module information
!-----------------------------------------------------------------------

use fms_mod,                       only : write_version_number

use zetac_horiz_grid_type_mod,     only : horiz_grid_type
use zetac_moisture_mod,            only : get_qsat
use zetac_vert_metric_type_mod,    only : vert_metric_type

!-----------------------------------------------------------------------
! public interfaces
!-----------------------------------------------------------------------

implicit none
private
public get_spec_hum

interface get_spec_hum
  module procedure get_spec_hum_yz
  module procedure get_spec_hum_xyz
end interface

character(len=*), parameter :: module='zetac_spec_hum_mod'
logical :: do_init=.true.

!-----------------------------------------------------------------------
! version control information
!-----------------------------------------------------------------------

character(len=128)  :: version =  '$Id: zetac_spec_hum.f90,v 1.1.2.3.2.1 2004/11/04 22:46:16 stg Exp $'
character(len=128)  :: tag     =  '$Name:  $'

contains

!#######################################################################

subroutine get_spec_hum_yz (                                           &
                                                       Mgrid, Vmetric, &
                                           ppref, taref, rhref, qvref )

type (horiz_grid_type),                        intent (in)    ::       &
                                                                Mgrid

type (vert_metric_type),                      intent (in)    ::        &
                                                              Vmetric

real, dimension(Mgrid%jbg:Mgrid%jeg,Mgrid%iz), intent (in)    ::       &
                                                  ppref, taref, rhref

real, dimension(Mgrid%jbg:Mgrid%jeg,Mgrid%iz), intent (out)   ::       &
                                                                qvref

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

integer :: j, jbg, jeg
integer :: k, kbd, ked


  jbg = Mgrid%jbg
  jeg = Mgrid%jeg
  kbd = Mgrid%kbd
  ked = Mgrid%ked

!-----------------------------------------------------------------------
! initialize module
!-----------------------------------------------------------------------

  if (do_init) then
     call spec_hum_init
  endif

!-----------------------------------------------------------------------
! get reference specific humidity from reference rh
!-----------------------------------------------------------------------

  call get_qsat ( ppref, taref, qvref )

  do k=kbd,ked
     do j=jbg,jeg
        qvref(j,k) = qvref(j,k)*rhref(j,k)
     enddo
  enddo

  return
end subroutine get_spec_hum_yz

!#######################################################################

subroutine get_spec_hum_xyz ( Mgrid, Vmetric, pp, ta, rh, qv )

type (horiz_grid_type),                        intent (in)    ::       &
                                                                Mgrid

type (vert_metric_type),                      intent (in)    ::        &
                                                              Vmetric

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),          intent (in)    ::       &
                                                           pp, ta, rh

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),          intent (out)   ::       &
                                                                   qv

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

integer :: i, ibd, ied
integer :: j, jbd, jed
integer :: k, kbd, ked

  ibd = Mgrid%ibd
  ied = Mgrid%ied
  jbd = Mgrid%jbd
  jed = Mgrid%jed
  kbd = Mgrid%kbd
  ked = Mgrid%ked

!-----------------------------------------------------------------------
! initialize module
!-----------------------------------------------------------------------

  if (do_init) then
     call spec_hum_init
  endif

!-----------------------------------------------------------------------
! get total specific humidity from total rh
!-----------------------------------------------------------------------

  call get_qsat ( pp, ta, qv )

  do k=kbd,ked
     do j=jbd,jed
        do i=ibd,ied
           qv(i,j,k) = qv(i,j,k)*rh(i,j,k)
        enddo
     enddo
  enddo

  return
end subroutine get_spec_hum_xyz

!#######################################################################

subroutine spec_hum_init

!-----------------------------------------------------------------------
! write version number to logfile
!-----------------------------------------------------------------------

  call write_version_number (version, tag)

  do_init = .false.
  
  return
end subroutine spec_hum_init

!#######################################################################

end module zetac_spec_hum_mod

