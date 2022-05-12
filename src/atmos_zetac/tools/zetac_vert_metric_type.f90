module zetac_vert_metric_type_mod

!-----------------------------------------------------------------------
! module information
!-----------------------------------------------------------------------

use fms_mod,                     only : write_version_number

use zetac_horiz_grid_type_mod,   only : horiz_grid_type

!-----------------------------------------------------------------------
! public interfaces
!-----------------------------------------------------------------------

implicit none
private
public vert_metric_type, vert_metric_type_define

character(len=*), parameter :: module='zetac_vert_metric_type_mod'

type vert_metric_type

real :: ptop, pbot

real, pointer, dimension(:)   :: zetaw   =>NULL()
real, pointer, dimension(:)   :: zetam   =>NULL()
real, pointer, dimension(:)   :: dzetaw  =>NULL()
real, pointer, dimension(:)   :: dzetam  =>NULL()
real, pointer, dimension(:,:) :: topog   =>NULL()

end type vert_metric_type

!-----------------------------------------------------------------------
! version control information
!-----------------------------------------------------------------------

character(len=128)  :: version =  '$Id: zetac_vert_metric_type.f90,v 1.1.2.3.2.1 2004/07/09 00:50:02 stg Exp $'
character(len=128)  :: tag     =  '$Name: latest $'

contains

!#######################################################################

subroutine vert_metric_type_define (                                   &
                                                       Mgrid, Vmetric, &
                                                           pbot, ptop, &
                                         zetaw, zetam, dzetam, dzetaw, &
                                                                topog )

type (horiz_grid_type),                        intent (in)    ::       &
                                                                Mgrid

type (vert_metric_type),                       intent (out)   ::       &
                                                              Vmetric

real   ,                                       intent (in)    ::       &
                                                           pbot, ptop

real   , dimension(Mgrid%iz),                  intent (in)    ::       &
                                         zetaw, zetam, dzetam, dzetaw

real   , dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed),          &
                                               intent (in)    ::       &
                                                                topog

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

  allocate (Vmetric%zetaw(Mgrid%iz))
  allocate (Vmetric%zetam(Mgrid%iz))
  allocate (Vmetric%dzetam(Mgrid%iz))
  allocate (Vmetric%dzetaw(Mgrid%iz))

  allocate (Vmetric%topog(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed))

!-----------------------------------------------------------------------
! write version number to logfile
!-----------------------------------------------------------------------

  call write_version_number (version, tag)

  Vmetric%pbot = pbot
  Vmetric%ptop = ptop

  Vmetric%zetaw  = zetaw
  Vmetric%zetam  = zetam
  Vmetric%dzetam = dzetam
  Vmetric%dzetaw = dzetaw

  Vmetric%topog = topog

  return
end subroutine vert_metric_type_define

!#######################################################################

end module zetac_vert_metric_type_mod
    

