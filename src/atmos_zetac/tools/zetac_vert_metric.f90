module zetac_vert_metric_mod

!-----------------------------------------------------------------------
! module information
!-----------------------------------------------------------------------

use fms_mod,                      only : write_version_number, stdout, &
                                         error_mesg, FATAL,mpp_pe

use zetac_horiz_grid_type_mod,    only : horiz_grid_type
use zetac_horiz_metric_type_mod,  only : horiz_metric_type
use zetac_vert_metric_type_mod,   only : vert_metric_type,             &
                                         vert_metric_type_define

!-----------------------------------------------------------------------
! public interfaces
!-----------------------------------------------------------------------

implicit none
private
public vert_metric

character(len=*), parameter :: module='zetac_vert_metric_mod'
logical :: do_init=.true.

!-----------------------------------------------------------------------
! version control information
!-----------------------------------------------------------------------

character(len=128)  :: version =  '$Id: zetac_vert_metric.f90,v 1.1.2.4.2.3 2004/07/09 23:16:21 stg Exp $'
character(len=128)  :: tag     =  '$Name: latest $'

contains

!#######################################################################

subroutine vert_metric (                                               &
                                              Mgrid, Hmetric, Vmetric, &
                                                           pbot, ptop, &
                                                         zetaw, topog )

type (horiz_grid_type),                     intent (in)    ::   Mgrid

type (horiz_metric_type),                   intent (in)    :: Hmetric

type (vert_metric_type),                    intent (out)   :: Vmetric

real,                                       intent (in)    ::          &
                                                           pbot, ptop

real, dimension(Mgrid%iz),                  intent (in)    ::          &
                                                                zetaw

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed),             &
                                            intent (in)    ::          &
                                                                topog

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

real, dimension(Mgrid%iz)                                   ::         &
                                                zetam, dzetam, dzetaw

integer :: k, iz

  iz = Mgrid%iz

!-----------------------------------------------------------------------
! initialize module
!-----------------------------------------------------------------------

  if (do_init) then
     call vert_metric_init
  endif

!-----------------------------------------------------------------------
! validate vertical grid
!-----------------------------------------------------------------------

  if (zetaw(1) /= 0. .or. zetaw(iz) /= 1.) then
print *,iz,zetaw
     call error_mesg (module, "invalid vertical grid", FATAL)
  endif

  do k=2,iz
     zetam(k) = 0.5*(zetaw(k) + zetaw(k-1))
  enddo
  zetam(1) = 0.

!-----------------------------------------------------------------------
! vertical grid intervals
!-----------------------------------------------------------------------

  do k=2,iz
     dzetam(k) = zetaw(k) - zetaw(k-1)
  enddo
  do k=2,iz-1
     dzetaw(k) = zetam(k+1) - zetam(k)
  enddo

  dzetam(1) = dzetam(2)
  dzetaw(1) = dzetam(2)
!  dzetaw(iz) = dzetam(iz)
  dzetaw(iz) = zetaw(iz) - zetam(iz)  !stg

  if (any(dzetam(2:iz) <= 0.)) then
     call error_mesg (module, "vertical grid is non-monotonic", FATAL)
  endif

!-----------------------------------------------------------------------
! define vertical metric derived type
!-----------------------------------------------------------------------
 
  call vert_metric_type_define (                                       &
                                                       Mgrid, Vmetric, &
                                                           pbot, ptop, &
                                         zetaw, zetam, dzetam, dzetaw, &
                                                                topog )

  write (stdout(), '(/,"zeta =",10f9.1/)') zetaw

  return
end subroutine vert_metric

!#######################################################################

subroutine vert_metric_init

!-----------------------------------------------------------------------
! write version number to logfile
!-----------------------------------------------------------------------

  call write_version_number (version, tag)

  do_init = .false.
  
  return
end subroutine vert_metric_init

!#######################################################################

end module zetac_vert_metric_mod
