module zetac_press_grad_mod

!-----------------------------------------------------------------------
! module information
!-----------------------------------------------------------------------

use mpp_mod,                       only : mpp_pe, mpp_root_pe
use fms_mod,                       only : write_version_number

use zetac_axes_mod,                only : gridv
use zetac_horiz_grid_type_mod,     only : horiz_grid_type
use zetac_horiz_metric_type_mod,   only : horiz_metric_type
use zetac_vert_metric_type_mod,    only : vert_metric_type

!-----------------------------------------------------------------------
! public interfaces
!-----------------------------------------------------------------------

implicit none
private
public get_pgf, press_grad_xy, press_grad_init

character(len=*), parameter :: module='zetac_press_grad_mod'

real, allocatable, dimension(:) :: dxu, zetaw

real :: dy
integer :: ibc, iec
integer :: jbc, jec
integer :: kbd, ked

!-----------------------------------------------------------------------
! version control information
!-----------------------------------------------------------------------

character(len=128)  :: version =  '$Id: zetac_press_grad.f90,v 1.1.2.3.2.3 2004/11/04 22:49:57 stg Exp $'
character(len=128)  :: tag     =  '$Name:  $'

contains

!#######################################################################

subroutine get_pgf ( Mgrid, gradp, alo, ahi, pgf )

type (horiz_grid_type),                     intent (in)       ::       &
                                                                Mgrid

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),       intent (in)       ::       &
                                                      gradp, alo, ahi

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),       intent (out)      ::  pgf

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

integer :: i, j, k

!-----------------------------------------------------------------------
! PGF = dpdx or dpdy times alpha
!-----------------------------------------------------------------------

  do k=kbd+1,ked
     do j=jbc-1,jec
        do i=ibc-1,iec
           pgf(i,j,k) = gradp(i,j,k-1)*ahi(i,j,k)                      &
                      + gradp(i,j,k)  *alo(i,j,k)
        enddo
     enddo
  enddo

  return
end subroutine get_pgf

!#######################################################################

subroutine press_grad_xy ( Mgrid, pp, dpdx, dpdy, fz, dfdx, dfdy )

type (horiz_grid_type),                         intent (in)  ::        &
                                                                Mgrid

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),           intent (in)   ::       &
                                                                   pp

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),           intent (out)  ::       &
                                                           dpdx, dpdy

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                          Mgrid%kbd:Mgrid%ked), intent (in)   ::       &
                                                                   fz

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                          Mgrid%kbd:Mgrid%ked), intent (out)   ::      &
                                                           dfdx, dfdy

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

integer :: i, j, k

!-----------------------------------------------------------------------
! horizontal pressure gradient
!-----------------------------------------------------------------------

  do k=kbd,ked
     do j=jbc-1,jec
        do i=ibc-1,iec
           dpdx(i,j,k) = (pp(i+1,j,k) - pp(i,j,k))/dxu(j)
           dpdy(i,j,k) = (pp(i,j+1,k) - pp(i,j,k))/dy
        enddo
     enddo
  enddo

  do k=kbd+1,ked
     do j=jbc-1,jec
        do i=ibc-1,iec
           dfdx(i,j,k) = (fz(i+1,j,k) - fz(i,j,k))/dxu(j)
           dfdy(i,j,k) = (fz(i,j+1,k) - fz(i,j,k))/dy
        enddo
     enddo
  enddo

  return
end subroutine press_grad_xy

!#######################################################################

subroutine press_grad_init ( Mgrid, Hmetric, Vmetric )

type (horiz_grid_type),                                                &
                                         intent (in)    ::      Mgrid

type (horiz_metric_type),                                              &
                                         intent (in)    ::    Hmetric

type (vert_metric_type),                 intent (in)    ::             &
                                                              Vmetric

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

  ibc  = Mgrid%ibc
  iec  = Mgrid%iec
  jbc  = Mgrid%jbc
  jec  = Mgrid%jec
  kbd  = Mgrid%kbd
  ked  = Mgrid%ked

!-----------------------------------------------------------------------
! allocate and define static arrays
!-----------------------------------------------------------------------

  allocate (dxu(Mgrid%jbd:Mgrid%jed))

  dxu    = Hmetric%dxu(Mgrid%jbd:Mgrid%jed)
  dy     = Hmetric%dy

  allocate (zetaw(Mgrid%iz))

  zetaw = Vmetric%zetaw

  return
end subroutine press_grad_init

!#######################################################################

end module zetac_press_grad_mod
