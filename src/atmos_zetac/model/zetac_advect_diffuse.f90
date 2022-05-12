module zetac_advect_diffuse_mod

!-----------------------------------------------------------------------
! module information
!-----------------------------------------------------------------------

use mpp_mod,                     only : mpp_pe, mpp_root_pe
use fms_mod,                     only : stdlog, write_version_number
use field_manager_mod,           only : model_atmos
use tracer_manager_mod,          only : get_tracer_index

use zetac_advect_horiz_mod,      only : advect_horiz_at_tracer
use zetac_advect_vert_mod,       only : advect_vert_at_tracer
use zetac_horiz_grid_type_mod,   only : horiz_grid_type
use zetac_horiz_metric_type_mod, only : horiz_metric_type
use zetac_mix_horiz_mod,         only : mix_horiz
use zetac_time_pointers_mod,     only : ntime
use zetac_tracer_mod,            only : tracer_type
use zetac_vert_metric_type_mod,  only : vert_metric_type

!-----------------------------------------------------------------------
! public interfaces
!-----------------------------------------------------------------------

implicit none
private
public advect_diffuse

character(len=*), parameter :: module='zetac_advect_diffuse_mod'
logical :: do_init=.true.

real, allocatable, dimension(:,:,:) :: advect_horiz, advect_vert,      &
                                       mix_tend

integer :: ibc, iec
integer :: jbc, jec
integer :: kbd, ked

!-----------------------------------------------------------------------
! version control information
!-----------------------------------------------------------------------

character(len=128)  :: version =  '$Id: zetac_advect_diffuse.f90,v 1.1.2.4.2.5 2005/07/05 02:46:43 stg Exp $'
character(len=128)  :: tag     =  '$Name:  $'

contains

!#######################################################################

subroutine advect_diffuse (                                            &
                                              Mgrid, Hmetric, Vmetric, & 
                                   uu, vv, oo, qc, dpu, dpv, dpm, dhm, &
                                                              qc_tend, &
                                         index, past, present, future, &
                                                                 delt )

type (horiz_grid_type),                         intent (in)    ::      &
                                                                Mgrid

type (horiz_metric_type),                       intent (in)    ::      &
                                                              Hmetric

type (vert_metric_type),                        intent (in)    ::      &
                                                              Vmetric

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked, ntime),    intent (in)    ::      &
                                        uu, vv, oo, qc, dpu, dpv, dpm

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),           intent (in)    ::      &
                                                                  dhm

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),           intent (out)   ::      &
                                                              qc_tend

real,                                           intent (in)    ::      &
                                                                 delt

integer,                                        intent (inout) ::      &
                                                                index

integer,                                        intent (in)    ::      &
                                                past, present, future

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

integer :: n
integer :: i, j, k

!-----------------------------------------------------------------------
! initialize module
!-----------------------------------------------------------------------

  if (do_init) then
     call advect_diffuse_init ( Mgrid )
  endif                                                

!-----------------------------------------------------------------------
! horizontal advection
!-----------------------------------------------------------------------

   call advect_horiz_at_tracer (                                       &
                                              Mgrid, Hmetric, Vmetric, &
                                            qc, uu, vv, dpu, dpv, dpm, &
                                                                index, &
                                                         advect_horiz, &
                                                                 delt, &
                                                past, present, future )

!-----------------------------------------------------------------------
! vertical advection
!-----------------------------------------------------------------------

   call advect_vert_at_tracer (                                        &
                                                       Mgrid, Vmetric, &
                                                     qc, oo, dpm, dhm, &
                                                                index, &
                                                          advect_vert, &
                                                                 delt, &
                                                past, present, future )

!-----------------------------------------------------------------------
! horizontal diffusion
!-----------------------------------------------------------------------

  call mix_horiz (                                                     &
                                                                Mgrid, &
                                                       qc(:,:,:,past), &
                                                             mix_tend, &
                                                                index, &
                                                 dp=dpm(:,:,:,future) )

!-----------------------------------------------------------------------
! explicit tendency
!-----------------------------------------------------------------------

  do k=kbd+1,ked
     do j=jbc,jec
        do i=ibc,iec
           qc_tend(i,j,k) = mix_tend(i,j,k)                            &
                         - (advect_horiz(i,j,k) + advect_vert(i,j,k))
        enddo
     enddo
  enddo

  return
end subroutine advect_diffuse

!#######################################################################

subroutine advect_diffuse_init ( Mgrid )

type (horiz_grid_type),                         intent (in)    ::      &
                                                                Mgrid

!-----------------------------------------------------------------------
! write version number and namelist to logfile
!-----------------------------------------------------------------------

  call write_version_number (version, tag)

  allocate (advect_horiz(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,     &
                         Mgrid%kbd:Mgrid%ked))
  allocate (advect_vert (Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,     &
                         Mgrid%kbd:Mgrid%ked))
  allocate (mix_tend    (Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,     &
                         Mgrid%kbd:Mgrid%ked))

  advect_horiz = 0.0
  advect_vert  = 0.0
  mix_tend = 0.0

  ibc = Mgrid%ibc
  iec = Mgrid%iec
  jbc = Mgrid%jbc
  jec = Mgrid%jec
  kbd = Mgrid%kbd
  ked = Mgrid%ked

  do_init = .false.

  return
end subroutine advect_diffuse_init

!#######################################################################

end module zetac_advect_diffuse_mod
