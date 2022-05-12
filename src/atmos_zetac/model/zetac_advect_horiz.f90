module zetac_advect_horiz_mod

!-----------------------------------------------------------------------
! module information
!-----------------------------------------------------------------------

use mpp_mod,                     only : mpp_pe, mpp_root_pe
use fms_mod,                     only : write_version_number
use vert_advection_mod,          only : ADVECTIVE_FORM, FLUX_FORM,     &
                                        SECOND_CENTERED,               &
                                        FOURTH_CENTERED,               &
                                        FINITE_VOLUME_LINEAR
use zonal_advection_mod,         only : zonal_advection
use merid_advection_mod,         only : merid_advection

use zetac_horiz_grid_type_mod,   only : horiz_grid_type
use zetac_horiz_metric_type_mod, only : horiz_metric_type
use zetac_time_pointers_mod,     only : ntime
use zetac_vert_metric_type_mod,  only : vert_metric_type
use zetac_tracer_mod,            only : get_method

!-----------------------------------------------------------------------
! public interfaces
!-----------------------------------------------------------------------

implicit none
private
public advect_horiz_at_u, advect_horiz_at_v,                           &
       advect_horiz_at_mass, advect_horiz_at_tracer,                   &
       advect_horiz_quick

character(len=*), parameter :: module='zetac_advect_horiz_mod'
character(len=16)    :: name_scheme, name_form
integer              :: num_scheme, num_form

real, allocatable, dimension(:)     :: dxu, dxv, cosulat, cosvlat
real, allocatable, dimension(:,:,:) :: uu1, vv1
real, allocatable, dimension(:,:,:) :: dx1, dy1
real, allocatable, dimension(:,:,:) :: var1, uvar, vvar
real, allocatable, dimension(:,:,:) :: advect_zonal, advect_merid
real, allocatable, dimension(:,:,:) :: densu, densv, udensu, vdensv

real    :: dy
integer :: ibc, iec, ibd, ied
integer :: jbc, jec, jbd, jed
integer :: kbc, kec, kbd, ked
integer :: nbuf

logical :: do_init=.true.

!-----------------------------------------------------------------------
! version control information
!-----------------------------------------------------------------------

character(len=128)  :: version =  '$Id: zetac_advect_horiz.f90,v 1.1.2.4.2.5 2005/06/16 18:00:37 stg Exp $'
character(len=128)  :: tag     =  '$Name:  $'

contains

!#######################################################################

subroutine advect_horiz_at_u (                                         &
                                              Mgrid, Hmetric, Vmetric, &
                                                  var, uu, vv, du, dv, &
                                                                index, &
                                                         advect_horiz, &
                                                        do_flux, delt, &
                                                           past, pres, &
                                                               advect )

!-----------------------------------------------------------------------
! calculate horizontal advection of var when var is at u-points
!-----------------------------------------------------------------------

type (horiz_grid_type),                                intent (in)  :: &
                                                                Mgrid

type (horiz_metric_type),                              intent (in)  :: &
                                                              Hmetric

type (vert_metric_type),                               intent (in)  :: &
                                                              Vmetric

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked, ntime),           intent (in)  :: &
                                                  var, uu, vv, du, dv

integer,                                               intent (in)  :: &
                                                                index

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),                  intent (out) :: &
                                                         advect_horiz

logical,                                               intent (in)  :: &
                                                              do_flux

real,                                                  intent (in)  :: &
                                                                 delt

integer,                                               intent (in)  :: &
                                                           past, pres

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked), optional,        intent (out) :: &
                                                               advect

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

real    :: dt
integer :: ltime
integer :: i, j, k

!-----------------------------------------------------------------------
! initialize module
!-----------------------------------------------------------------------

  if (do_init) then
     call advect_horiz_init ( Mgrid, Hmetric, Vmetric )
  endif

  call get_method ( index, 'advect_horiz', name_scheme, num_scheme )
  call get_method ( index, 'equation', name_form, num_form )

  if ( nbuf < 2 ) num_scheme = min(num_scheme, SECOND_CENTERED)

  if ( num_scheme == SECOND_CENTERED .and.                             &
       name_scheme(1:3) == "alt" .or.                                  &
       present(advect) ) then

!-----------------------------------------------------------------------
! second-centered advective form (inline)
!-----------------------------------------------------------------------

     densv = dv(:,:,:,pres)
     densu = du(:,:,:,pres)
     vdensv = vv(:,:,:,pres)*densv
     udensu = uu(:,:,:,pres)*densu

     if ( do_flux ) then

        do j=jbc,jec
           do i=ibc,iec+1
              uvar(i,j,:) = 0.25*(udensu(i-1,j,:) + udensu(i,j,:))*    &
                                 (var(i-1,j,:,pres) + var(i,j,:,pres))
           enddo
        enddo
        do j=jbc-1,jec
           do i=ibc,iec
              vvar(i,j,:) = 0.25*(vdensv(i,j,:) + vdensv(i+1,j,:))*    &
                                 (var(i,j,:,pres) + var(i,j+1,:,pres))
           enddo
        enddo

        do k=kbd+1,ked
           do j=jbc,jec
              do i=ibc,iec
                  advect_horiz(i,j,k) =                                &
                         ((uvar(i+1,j,k) - uvar(i,j,k))/dxu(j)         &
                        + (vvar(i,j,k) - vvar(i,j-1,k))/dy)/           &
                                                          densu(i,j,k)
    
              enddo
           enddo
        enddo

     else

       do k=kbd+1,ked
           do j=jbc,jec
              do i=ibc,iec
                 advect_horiz(i,j,k) = 0.25*(                          &
                   ( (uu(i+1,j,k,pres) +  uu(i,j,k,pres))*             &
                    (var(i+1,j,k,pres) - var(i,j,k,pres))              &
                     + (uu(i,j,k,pres) +  uu(i-1,j,k,pres))*           &
                      (var(i,j,k,pres) - var(i-1,j,k,pres)) ) / dxu(j) &
                 + ( (vv(i+1,j,k,pres) +  vv(i,j,k,pres))*             &
                    (var(i,j+1,k,pres) - var(i,j,k,pres))              &
                 + (vv(i+1,j-1,k,pres) +  vv(i,j-1,k,pres))*           &
                      (var(i,j,k,pres) - var(i,j-1,k,pres)) ) / dy )
              enddo
           enddo
        enddo

     endif

     if (present(advect)) advect = advect_horiz

     if ( num_scheme == SECOND_CENTERED .and.                          &
       name_scheme(1:3) == "alt" ) return

  endif

!-----------------------------------------------------------------------
! all other cases
!-----------------------------------------------------------------------

  if ( num_scheme <= FOURTH_CENTERED ) then
     ltime = pres
  else 
     ltime = past
  endif

  do k=kbd+1,ked  ! note: uu1/vv1 are dimensioned from ibd-1/jbd-1
     do j=jbd,jed
        do i=ibd,ied-1
           uu1(i,j,k) = 0.5*(uu(i,j,k,pres) + uu(i+1,j,k,pres))
           vv1(i,j,k) = 0.5*(vv(i,j,k,pres) + vv(i+1,j,k,pres))
        enddo
     enddo
  enddo

  do k=kbd+1,ked
     do j=jbd,jed
        uu1(ied,j,k) = uu(ied,j,k,pres)
        vv1(ied,j,k) = vv(ied,j,k,pres)
     enddo
  enddo

  do k=kbd+1,ked
     do j=jbd,jed
        uu1(ibd-1,j,k) = uu1(ibd,j,k)
     enddo
     do i=ibd,ied
        vv1(i,jbd-1,k) = vv1(i,jbd,k)
     enddo
  enddo

  do k=kbd+1,ked
     do j=jbd,jed
        do i=ibd,ied
           var1(i,j,k) = var(i,j,k,ltime)
            dx1(i,j,k) = dxu(j)
            dy1(i,j,k) = dy
        enddo
     enddo
  enddo

  if ( do_flux ) then

     do k=kbd+1,ked
        do j=jbd,jed
           do i=ibd,ied
              var1(i,j,k) = var1(i,j,k)*du(i,j,k,ltime)
           enddo
        enddo
     enddo

  endif

  dt = delt

!-----------------------------------------------------------------------
! zonal advection
!-----------------------------------------------------------------------

  if ( Mgrid%ix > 1 ) then
     call zonal_advection ( dt, uu1, dx1, var1, advect_zonal,          &
                            scheme=num_scheme, form=FLUX_FORM )
  else
     advect_zonal = 0.0
  endif

!-----------------------------------------------------------------------
! meridional advection
!-----------------------------------------------------------------------

  if ( Mgrid%iy > 1 ) then

     if ( num_form == FLUX_FORM ) then
       do k=kbd+1,ked
           do j=jbd,jed
              do i=ibd,ied
                 var1(i,j,k) = var1(i,j,k)*cosulat(j)
              enddo
           enddo
        enddo
     endif

     call merid_advection ( dt, vv1, dy1, var1, advect_merid,          &
                            scheme=num_scheme, form=FLUX_FORM )

  else
     advect_merid = 0.0
  endif

  if ( do_flux ) then

!-----------------------------------------------------------------------
! complete flux-form tendency
!-----------------------------------------------------------------------

     do k=kbd+1,ked
        do j=jbc,jec
           do i=ibc,iec
              advect_horiz(i,j,k) =                                    &
              -(advect_zonal(i,j,k) + advect_merid(i,j,k)/cosulat(j))/ &
                                                         du(i,j,k,pres)
           enddo
        enddo
     enddo

  else

!-----------------------------------------------------------------------
! change back to advect form
!-----------------------------------------------------------------------

     do k=kbd+1,ked
        do j=jbc,jec
           do i=ibc,iec
              advect_horiz(i,j,k) =                                    &
                          -(advect_zonal(i,j,k) + advect_merid(i,j,k)) &
                 - var1(i,j,k)*((uu1(i,j,k) - uu1(i-1,j,k))/dx1(i,j,k) &
                              + (vv1(i,j,k) - vv1(i,j-1,k))/dy1(i,j,k))
           enddo
        enddo
     enddo

  endif

  return
end subroutine advect_horiz_at_u

!#######################################################################

subroutine advect_horiz_at_v (                                         &
                                              Mgrid, Hmetric, Vmetric, &
                                                  var, uu, vv, du, dv, &
                                                                index, &
                                                         advect_horiz, &
                                                        do_flux, delt, &
                                                           past, pres, &
                                                               advect )

!-----------------------------------------------------------------------
! calculate horizontal advection of var when var is at v-points 
!-----------------------------------------------------------------------
       
type (horiz_grid_type),                                intent (in)  :: &
                                                                Mgrid

type (horiz_metric_type),                              intent (in)  :: &
                                                              Hmetric

type (vert_metric_type),                               intent (in)  :: &
                                                              Vmetric

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked, ntime),           intent (in)  :: &
                                                  var, uu, vv, du, dv

integer,                                               intent (in)  :: &
                                                                index

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),                  intent (out) :: &

                                                         advect_horiz

logical,                                               intent (in)  :: &
                                                              do_flux

real,                                                  intent (in)  :: &
                                                                 delt

integer,                                               intent (in)  :: &
                                                           past, pres

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked), optional,        intent (out) :: &
                                                               advect

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------
       
real    :: dt
integer :: ltime
integer :: i, j, k

!-----------------------------------------------------------------------
! initialize module
!-----------------------------------------------------------------------

  if (do_init) then
     call advect_horiz_init ( Mgrid, Hmetric, Vmetric )
  endif

  call get_method ( index, 'advect_horiz', name_scheme, num_scheme )
  call get_method ( index, 'equation', name_form, num_form )

  if ( nbuf < 2 ) num_scheme = min(num_scheme, SECOND_CENTERED)

  if ( num_scheme == SECOND_CENTERED .and.                             &
       name_scheme(1:3) == "alt" .or.                                  &
       present(advect) ) then

!-----------------------------------------------------------------------
! second-centered advective form (inline)
!-----------------------------------------------------------------------

     densv = dv(:,:,:,pres)
     densu = du(:,:,:,pres)
     vdensv = vv(:,:,:,pres)*densv
     udensu = uu(:,:,:,pres)*densu

     if ( do_flux ) then

        do j=jbc,jec+1
           do i=ibc,iec
              vvar(i,j,:) = 0.25*(vdensv(i,j-1,:) + vdensv(i,j,:))*   &
                                 (var(i,j-1,:,pres) + var(i,j,:,pres))
           enddo
        enddo
        do j=jbc,jec
           do i=ibc-1,iec
              uvar(i,j,:) = 0.25*(udensu(i,j,:) + udensu(i,j+1,:))*   &
                                 (var(i,j,:,pres) + var(i+1,j,:,pres))
           enddo
        enddo

        do k=kbd+1,ked
           do j=jbc,jec
              do i=ibc,iec
                  advect_horiz(i,j,k) =                                &
                          ((vvar(i,j+1,k) - vvar(i,j,k))/dy            &
                         + (uvar(i,j,k) - uvar(i-1,j,k))/dxv(j))/      &
                                                          densv(i,j,k)
    
              enddo
           enddo
        enddo

     else

        do k=kbd+1,ked
           do j=jbc,jec
              do i=ibc,iec
                 advect_horiz(i,j,k) = 0.25*(                          &
                   ( (uu(i,j+1,k,pres) +  uu(i,j,k,pres))*             &
                    (var(i+1,j,k,pres) - var(i,j,k,pres))              &
                 + (uu(i-1,j+1,k,pres) +  uu(i-1,j,k,pres))*           &
                      (var(i,j,k,pres) - var(i-1,j,k,pres)) ) / dxv(j) &
                 + ( (vv(i,j+1,k,pres) +  vv(i,j,k,pres))*             &
                    (var(i,j+1,k,pres) - var(i,j,k,pres))              &
                     + (vv(i,j,k,pres) +  vv(i,j-1,k,pres))*           &
                      (var(i,j,k,pres) - var(i,j-1,k,pres)) ) / dy )
              enddo
           enddo
        enddo

     endif

     if (present(advect)) advect = advect_horiz

     if ( num_scheme == SECOND_CENTERED .and.                          &
       name_scheme(1:3) == "alt" ) return

  endif

!-----------------------------------------------------------------------
! all other cases
!-----------------------------------------------------------------------

  if ( num_scheme <= FOURTH_CENTERED ) then
     ltime = pres
  else 
     ltime = past
  endif

  do k=kbd+1,ked  ! note: uu1/vv1 are dimensioned from ibd-1 & jbd-1
     do j=jbd,jed-1
        do i=ibd,ied
           uu1(i,j,k) = 0.5*(uu(i,j,k,pres) + uu(i,j+1,k,pres))
           vv1(i,j,k) = 0.5*(vv(i,j,k,pres) + vv(i,j+1,k,pres))
        enddo
     enddo
  enddo

  do k=kbd+1,ked
     do i=ibd,ied
        uu1(i,jed,k) = uu(i,jed,k,pres)
        vv1(i,jed,k) = vv(i,jed,k,pres)
     enddo
  enddo

  do k=kbd+1,ked
     do j=jbd,jed
        uu1(ibd-1,j,k) = uu1(ibd,j,k)
     enddo
     do i=ibd,ied
        vv1(i,jbd-1,k) = vv1(i,jbd,k)
     enddo
  enddo

  do k=kbd+1,ked
     do j=jbd,jed
        do i=ibd,ied
           var1(i,j,k) = var(i,j,k,ltime)
            dx1(i,j,k) = dxv(j)
            dy1(i,j,k) = dy
        enddo
     enddo
  enddo

  if ( do_flux ) then

     do k=kbd+1,ked
        do j=jbd,jed
           do i=ibd,ied
              var1(i,j,k) = var1(i,j,k)*dv(i,j,k,ltime)
           enddo
        enddo
     enddo

  endif

  dt = delt

!-----------------------------------------------------------------------
! zonal advection
!-----------------------------------------------------------------------

  if ( Mgrid%ix > 1 ) then
     call zonal_advection ( dt, uu1, dx1, var1, advect_zonal,          &
                            scheme=num_scheme, form=FLUX_FORM )
  else
     advect_zonal = 0.0
  endif

!-----------------------------------------------------------------------
! meridional advection
!-----------------------------------------------------------------------

  if ( Mgrid%iy > 1 ) then

     if ( num_form == FLUX_FORM ) then
        do k=kbd+1,ked
           do j=jbd,jed
              do i=ibd,ied
                 var1(i,j,k) = var1(i,j,k)*cosvlat(j)
              enddo
           enddo
        enddo
     endif

     call merid_advection ( dt, vv1, dy1, var1, advect_merid,          &
                            scheme=num_scheme, form=FLUX_FORM )

  else
     advect_merid = 0.0
  endif

  if ( do_flux ) then

!-----------------------------------------------------------------------
! complete flux-form tendency
!-----------------------------------------------------------------------

     do k=kbd+1,ked
        do j=jbc,jec
           do i=ibc,iec
              advect_horiz(i,j,k) =                                    &
              -(advect_zonal(i,j,k) + advect_merid(i,j,k)/cosvlat(j))/ &
                                                        dv(i,j,k,pres)
           enddo
        enddo
     enddo

  else

!-----------------------------------------------------------------------
! change back to advect form
!-----------------------------------------------------------------------

     do k=kbd+1,ked
        do j=jbc,jec
           do i=ibc,iec
              advect_horiz(i,j,k) =                                    &
                          -(advect_zonal(i,j,k) + advect_merid(i,j,k)) &
                 - var1(i,j,k)*((uu1(i,j,k) - uu1(i-1,j,k))/dx1(i,j,k) &
                              + (vv1(i,j,k) - vv1(i,j-1,k))/dy1(i,j,k))
           enddo
        enddo
     enddo

  endif

  return
end subroutine advect_horiz_at_v

!#######################################################################

subroutine advect_horiz_at_mass (                                      &
                                              Mgrid, Hmetric, Vmetric, &
                                              var, uu, vv, du, dv, dm, &
                                                                index, &
                                                         advect_horiz, &
                                                        do_flux, delt, &
                                                           past, pres, &
                                                               advect )

!-----------------------------------------------------------------------
! calculate horizontal advection of var when var is at mass-points 
!-----------------------------------------------------------------------

type (horiz_grid_type),                                intent (in)  :: &
                                                                Mgrid

type (horiz_metric_type),                              intent (in)  :: &
                                                              Hmetric

type (vert_metric_type),                               intent (in)  :: &
                                                              Vmetric

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked, ntime),           intent (in)  :: &
                                              var, uu, vv, du, dv, dm

integer,                                               intent (in)  :: &
                                                                index

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),                  intent (out) :: &
                                                         advect_horiz

logical,                                               intent (in)  :: &
                                                              do_flux

real,                                                  intent (in)  :: &
                                                                 delt

integer,                                               intent (in)  :: &
                                                           past, pres

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked), optional,        intent (out) :: &
                                                               advect

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

real    :: dt
integer :: ltime
integer :: i, j, k

!-----------------------------------------------------------------------
! initialize module
!-----------------------------------------------------------------------

  if (do_init) then
     call advect_horiz_init ( Mgrid, Hmetric, Vmetric )
  endif

  call get_method ( index, 'advect_horiz', name_scheme, num_scheme )
  call get_method ( index, 'equation', name_form, num_form )

  if ( nbuf < 2 ) num_scheme = min(num_scheme, SECOND_CENTERED)

  if ( num_scheme == SECOND_CENTERED .and.                             &
       name_scheme(1:3) == "alt" .or.                                  &
       present(advect) ) then

!-----------------------------------------------------------------------
! second-centered (inline)
!-----------------------------------------------------------------------

     densv = dv(:,:,:,pres)
     densu = du(:,:,:,pres)
     vdensv = vv(:,:,:,pres)*densv
     udensu = uu(:,:,:,pres)*densu

     if ( do_flux ) then

        do j=jbc-1,jec
           do i=ibc,iec
              vvar(i,j,:) = 0.5*vdensv(i,j,:)*                         &
                                 (var(i,j,:,pres) + var(i,j+1,:,pres))
           enddo
        enddo
        do j=jbc,jec
           do i=ibc-1,iec
              uvar(i,j,:) = 0.5*udensu(i,j,:)*                         &
                                 (var(i,j,:,pres) + var(i+1,j,:,pres))
           enddo
        enddo

        do k=kbd+1,ked
           do j=jbc,jec
              do i=ibc,iec
                  advect_horiz(i,j,k) = (                              &
                           (uvar(i,j,k) - uvar(i-1,j,k))/dxu(j)        &
                         + (vvar(i,j,k) - vvar(i,j-1,k))/dy            &
                                                    ) / dm(i,j,k,pres)   
              enddo
           enddo
        enddo

     else

        do k=kbd+1,ked
           do j=jbc,jec
              do i=ibc,iec
                 advect_horiz(i,j,k) = 0.5*(                            &
                   ( udensu(i,j,k)*                                     &
                 (var(i+1,j,k,pres) - var(i,j,k,pres))                  &
                 + udensu(i-1,j,k)*                                     &
                   (var(i,j,k,pres) - var(i-1,j,k,pres)) ) / dxu(j)     &
                 + ( vdensv(i,j,k)*                                     &
                 (var(i,j+1,k,pres) - var(i,j,k,pres))                  &
                 + vdensv(i,j-1,k)*                                     &
                 (var(i,j,k,pres) - var(i,j-1,k,pres)) ) / dy           &
                                                     ) / dm(i,j,k,pres)
              enddo
           enddo
        enddo

     endif

     if (present(advect)) advect = advect_horiz

     if ( num_scheme == SECOND_CENTERED .and.                          &
       name_scheme(1:3) == "alt" ) return

  endif

!-----------------------------------------------------------------------
! all other cases
!-----------------------------------------------------------------------

  if ( num_scheme <= FOURTH_CENTERED ) then
     ltime = pres
  else 
     ltime = past
  endif

  do k=kbd+1,ked
     do j=jbd,jed
        do i=ibd,ied
           uu1(i,j,k) = uu(i,j,k,pres)
           vv1(i,j,k) = vv(i,j,k,pres)
        enddo
     enddo
  enddo

  do k=kbd+1,ked
     do j=jbd,jed
        uu1(ibd-1,j,k) = uu1(ibd,j,k)
     enddo
     do i=ibd,ied
        vv1(i,jbd-1,k) = vv1(i,jbd,k)
     enddo
  enddo

  do k=kbd+1,ked
     do j=jbd,jed
        do i=ibd,ied
           var1(i,j,k) = var(i,j,k,ltime)
            dx1(i,j,k) = dxu(j)
            dy1(i,j,k) = dy
        enddo
     enddo
  enddo

  if ( do_flux ) then

     do k=kbd+1,ked
        do j=jbd,jed
           do i=ibd,ied
              var1(i,j,k) = var1(i,j,k)*dm(i,j,k,ltime)
           enddo
        enddo
     enddo
     
  endif

  dt = delt

!-----------------------------------------------------------------------
! zonal advection
!-----------------------------------------------------------------------

  if ( Mgrid%ix > 1 ) then
     call zonal_advection ( dt, uu1, dx1, var1, advect_zonal,          &
                            scheme=num_scheme, form=FLUX_FORM )
  else
     advect_zonal = 0.0
  endif

!-----------------------------------------------------------------------
! meridional advection
!-----------------------------------------------------------------------

  if ( Mgrid%iy > 1 ) then

     if ( do_flux ) then
        do k=kbd+1,ked
           do j=jbd,jed
              do i=ibd,ied
                 var1(i,j,k) = var1(i,j,k)*cosulat(j)
              enddo
           enddo
        enddo
     endif

     call merid_advection ( dt, vv1, dy1, var1, advect_merid,          &
                            scheme=num_scheme, form=FLUX_FORM )

  else
     advect_merid = 0.0
  endif

  if ( do_flux ) then

!-----------------------------------------------------------------------
! complete flux-form tendency
!-----------------------------------------------------------------------

     do k=kbd+1,ked
        do j=jbc,jec
           do i=ibc,iec
              advect_horiz(i,j,k) =                                    &
              -(advect_zonal(i,j,k) + advect_merid(i,j,k)/cosulat(j))/ &
                                                         dm(i,j,k,pres)
           enddo
        enddo
     enddo

  else

!-----------------------------------------------------------------------
! change back to advect form
!-----------------------------------------------------------------------

     do k=kbd+1,ked
        do j=jbc,jec
           do i=ibc,iec
              advect_horiz(i,j,k) =                                    &
                          -(advect_zonal(i,j,k) + advect_merid(i,j,k)) &
                 - var1(i,j,k)*((uu1(i,j,k) - uu1(i-1,j,k))/dx1(i,j,k) &
                              + (vv1(i,j,k) - vv1(i,j-1,k))/dy1(i,j,k))
           enddo
        enddo
     enddo

  endif

  return
end subroutine advect_horiz_at_mass

!#######################################################################

subroutine advect_horiz_at_tracer (                                    &
                                              Mgrid, Hmetric, Vmetric, &
                                              var, uu, vv, du, dv, dm, &
                                                                index, &
                                                         advect_horiz, &
                                                                 delt, &
                                                   past, pres, future )

type (horiz_grid_type),                         intent (in)    ::      &
                                                                Mgrid

type (horiz_metric_type),                       intent (in)    ::      &
                                                              Hmetric

type (vert_metric_type),                               intent (in)  :: &
                                                              Vmetric

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked, ntime),    intent (in)    ::      &
                                              var, uu, vv, du, dv, dm

integer,                                        intent (inout) ::       &
                                                                index

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),           intent (out)   ::      &
                                                         advect_horiz

real,                                           intent (in)    ::      &
                                                                 delt

integer,                                        intent (in)    ::      &
                                                   past, pres, future

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

real    :: dt
integer :: ltime
integer :: i, j, k

!-----------------------------------------------------------------------
! initialize module
!-----------------------------------------------------------------------

  if (do_init) then
     call advect_horiz_init ( Mgrid, Hmetric, Vmetric )
  endif
  
  call get_method ( index, 'advect_horiz', name_scheme, num_scheme )
  call get_method ( index, 'equation', name_form, num_form )

  if ( name_scheme(1:4) == 'none' ) then
     do k=kbd+1,ked
        do j=jbc,jec
           do i=ibc,iec
              advect_horiz(i,j,k) = 0.0
           enddo
        enddo
     enddo
     return
  endif

  if ( nbuf < 2 ) num_scheme = min(num_scheme, SECOND_CENTERED)

  if ( num_scheme == SECOND_CENTERED .and.                             &
       name_scheme(1:3) == "alt" ) then

!-----------------------------------------------------------------------
! second-centered advective form (inline)
!-----------------------------------------------------------------------

     densv = dv(:,:,:,pres)
     densu = du(:,:,:,pres)
     vdensv = vv(:,:,:,pres)*densv
     udensu = uu(:,:,:,pres)*densu

     if ( num_form == FLUX_FORM ) then

        do j=jbc-1,jec
           do i=ibc,iec
              vvar(i,j,:) = 0.5*vdensv(i,j,:)*                         &
                                 (var(i,j,:,pres) + var(i,j+1,:,pres))
           enddo
        enddo
        do j=jbc,jec
           do i=ibc-1,iec
              uvar(i,j,:) = 0.5*udensu(i,j,:)*                         &
                                 (var(i,j,:,pres) + var(i+1,j,:,pres))
           enddo
        enddo

        do k=kbd+1,ked
           do j=jbc,jec
              do i=ibc,iec
                  advect_horiz(i,j,k) =                                &
                          ((uvar(i,j,k) - uvar(i-1,j,k))/dxv(j)        &
                         + (vvar(i,j,k) - vvar(i,j-1,k))/dy)/          &
                                                        dm(i,j,k,pres)
    
              enddo
           enddo
        enddo

     else

        do k=kbd+1,ked
           do j=jbc,jec
              do i=ibc,iec
                 advect_horiz(i,j,k) = 0.5*(                           &
                   ( uu(i,j,k,pres)*                                   &
                 (var(i+1,j,k,pres) - var(i,j,k,pres))                 &
                 + uu(i-1,j,k,pres)*                                   &
                   (var(i,j,k,pres) - var(i-1,j,k,pres)) ) / dxu(j)    &
                 + ( vv(i,j,k,pres)*                                   &
                 (var(i,j+1,k,pres) - var(i,j,k,pres))                 &
                 + vv(i,j-1,k,pres)*                                   &
                   (var(i,j,k,pres) - var(i,j-1,k,pres)) ) / dy )
              enddo
           enddo
        enddo

     endif

     return

  endif

!-----------------------------------------------------------------------
! all other cases
!-----------------------------------------------------------------------

  if ( num_scheme <= FOURTH_CENTERED ) then
     ltime = pres
  else
     ltime = past
  endif

  do k=kbd+1,ked
     do j=jbd,jed
        do i=ibd,ied
           uu1(i,j,k) = uu(i,j,k,pres)
           vv1(i,j,k) = vv(i,j,k,pres)
        enddo
     enddo
  enddo

  do k=kbd+1,ked
     do j=jbd,jed
        uu1(ibd-1,j,k) = uu1(ibd,j,k)
     enddo
     do i=ibd,ied
        vv1(i,jbd-1,k) = vv1(i,jbd,k)
     enddo
  enddo

  do k=kbd+1,ked
     do j=jbd,jed
        do i=ibd,ied
           var1(i,j,k) = var(i,j,k,ltime)
            dx1(i,j,k) = dxu(j)
            dy1(i,j,k) = dy
        enddo
     enddo
  enddo

  if ( num_form == FLUX_FORM ) then

     do k=kbd+1,ked
        do j=jbd,jed
           do i=ibd,ied
              var1(i,j,k) = var1(i,j,k)*dm(i,j,k,ltime)
           enddo
        enddo
     enddo

  endif

  dt = delt

!-----------------------------------------------------------------------
! zonal advection
!-----------------------------------------------------------------------

  if ( Mgrid%ix > 1 ) then
     call zonal_advection ( dt, uu1, dx1, var1, advect_zonal,          &
                            scheme=num_scheme, form=FLUX_FORM )
  else
     advect_zonal = 0.0
  endif

!-----------------------------------------------------------------------
! meridional advection
!-----------------------------------------------------------------------

  if ( Mgrid%iy > 1 ) then

     if ( num_form == FLUX_FORM ) then
        do k=kbd+1,ked
           do j=jbd,jed
              do i=ibd,ied
                 var1(i,j,k) = var1(i,j,k)*cosulat(j)
              enddo
           enddo
        enddo
     endif

     call merid_advection ( dt, vv1, dy1, var1, advect_merid,          &
                            scheme=num_scheme, form=FLUX_FORM )

  else
     advect_merid = 0.0
  endif

  if ( num_form == FLUX_FORM ) then

!-----------------------------------------------------------------------
! complete flux-form tendency
!-----------------------------------------------------------------------

     do k=kbd+1,ked
        do j=jbc,jec
           do i=ibc,iec
              advect_horiz(i,j,k) =                                    &
              -(advect_zonal(i,j,k) + advect_merid(i,j,k)/cosulat(j))/ &
                                                       dm(i,j,k,future)
           enddo
        enddo
     enddo

  else

!-----------------------------------------------------------------------
! change back to advect form
!-----------------------------------------------------------------------

     do k=kbd+1,ked
        do j=jbc,jec
           do i=ibc,iec
              advect_horiz(i,j,k) =                                    &
                          -(advect_zonal(i,j,k) + advect_merid(i,j,k)) &
                 - var1(i,j,k)*((uu1(i,j,k) - uu1(i-1,j,k))/dx1(i,j,k) &
                              + (vv1(i,j,k) - vv1(i,j-1,k))/dy1(i,j,k))
           enddo
        enddo
     enddo

  endif

  return
end subroutine advect_horiz_at_tracer

!#######################################################################

subroutine advect_horiz_quick (                                        &
                                              Mgrid, Hmetric, Vmetric, &
                                              var, uu, vv, du, dv, dm, &
                                                         advect_horiz )

!-----------------------------------------------------------------------
! calculate horizontal advection of var when var is at mass-points 
!-----------------------------------------------------------------------

type (horiz_grid_type),                                intent (in)  :: &
                                                                Mgrid

type (horiz_metric_type),                              intent (in)  :: &
                                                              Hmetric

type (vert_metric_type),                               intent (in)  :: &
                                                              Vmetric

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),                  intent (in)  :: &
                                              var, uu, vv, du, dv, dm

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),                  intent (out) :: &
                                                         advect_horiz

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

integer :: i, j, k

!-----------------------------------------------------------------------
! initialize module
!-----------------------------------------------------------------------

  if (do_init) then
     call advect_horiz_init ( Mgrid, Hmetric, Vmetric )
  endif

!-----------------------------------------------------------------------
! second-centered
!-----------------------------------------------------------------------

  vdensv = vv*dv
  udensu = uu*du

  do k=kbd+1,ked
     do j=jbc,jec
        do i=ibc,iec
           advect_horiz(i,j,k) =                                       &
            .5*(( udensu(i,j,k)*(var(i+1,j,k) - var(i,j,k))            &
              + udensu(i-1,j,k)*(var(i,j,k) - var(i-1,j,k)) ) / dxu(j) &
              + ( vdensv(i,j,k)*(var(i,j+1,k) - var(i,j,k))            &
              + vdensv(i,j-1,k)*(var(i,j,k) - var(i,j-1,k)) ) / dy )   &
                                                            / dm(i,j,k)
        enddo
     enddo
  enddo

end subroutine advect_horiz_quick

!#######################################################################

subroutine advect_horiz_init ( Mgrid, Hmetric, Vmetric )

type (horiz_grid_type),                             intent (in)  ::    &
                                                                Mgrid

type (horiz_metric_type),                           intent (in)  ::    &
                                                              Hmetric

type (vert_metric_type),                               intent (in)  :: &
                                                              Vmetric

  ibc = Mgrid%ibc
  iec = Mgrid%iec
  jbc = Mgrid%jbc
  jec = Mgrid%jec

  ibd = Mgrid%ibd
  ied = Mgrid%ied
  jbd = Mgrid%jbd
  jed = Mgrid%jed

  kbc = Mgrid%kbc
  kec = Mgrid%kec
  kbd = Mgrid%kbd
  ked = Mgrid%ked

  nbuf = Mgrid%nbuf

!-----------------------------------------------------------------------
! write version number to logfile
!-----------------------------------------------------------------------

  call write_version_number (version, tag)

!-----------------------------------------------------------------------
! allocate work arrays
!-----------------------------------------------------------------------

  allocate (uu1(ibd-1:ied,jbd:jed,kbc:kec))
  allocate (vv1(ibd:ied,jbd-1:jed,kbc:kec))

  allocate (var1(ibd:ied,jbd:jed,kbc:kec))
  allocate ( dx1(ibd:ied,jbd:jed,kbc:kec))
  allocate ( dy1(ibd:ied,jbd:jed,kbc:kec))

  allocate (advect_zonal(ibd:ied,jbd:jed,kbc:kec))
  allocate (advect_merid(ibd:ied,jbd:jed,kbc:kec))

  allocate (densu (ibd:ied,jbd:jed,kbd:ked))
  allocate (densv (ibd:ied,jbd:jed,kbd:ked))
  allocate (udensu(ibd:ied,jbd:jed,kbd:ked))
  allocate (vdensv(ibd:ied,jbd:jed,kbd:ked))
  allocate (uvar  (ibd:ied,jbd:jed,kbd:ked))
  allocate (vvar  (ibd:ied,jbd:jed,kbd:ked))

  allocate (dxu(jbd:jed))
  allocate (dxv(jbd:jed))
  allocate (cosulat(jbd:jed))
  allocate (cosvlat(jbd:jed))

  dy  = Hmetric%dy
  dxu = Hmetric%dxu(jbd:jed)
  dxv = Hmetric%dxv(jbd:jed)
  cosulat = Hmetric%cosulat(jbd:jed)
  cosvlat = Hmetric%cosvlat(jbd:jed)

  do_init = .false.

  return
end subroutine advect_horiz_init

!#######################################################################

end module zetac_advect_horiz_mod
