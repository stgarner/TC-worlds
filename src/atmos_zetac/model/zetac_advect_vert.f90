module zetac_advect_vert_mod

!-----------------------------------------------------------------------
! module information
!-----------------------------------------------------------------------

use mpp_mod,                     only : mpp_pe, mpp_root_pe
use fms_mod,                     only : write_version_number
use vert_advection_mod,          only : vert_advection,                &
                                        ADVECTIVE_FORM, FLUX_FORM,     &
                                        SECOND_CENTERED,               &
                                        FOURTH_CENTERED

use zetac_field_names_mod
use zetac_horiz_grid_type_mod,   only : horiz_grid_type
use zetac_time_pointers_mod,     only : ntime
use zetac_tracer_mod,            only : get_method
use zetac_vert_metric_type_mod,  only : vert_metric_type

!-----------------------------------------------------------------------
! public interfaces
!-----------------------------------------------------------------------

implicit none
private
public advect_vert_at_u, advect_vert_at_v,                             &
       advect_vert_at_mass, advect_vert_at_tracer,                     &
       advect_vert_quick

character(len=*), parameter :: module='zetac_advect_vert_mod'
character(len=16)  :: name_scheme, name_form
integer            :: num_scheme, num_form

integer :: ibc, iec
integer :: jbc, jec
integer :: kbd, ked
logical :: do_init=.true.

real, allocatable, dimension(:,:,:) :: ww_fms, dens_fms, fall_fms,     &
                                       var_fms, dz_fms, advect_fms,    &
                                       watuv, wvar

real, allocatable, dimension(:) :: dzetaw, dzetam

real, parameter :: addfac=0.4, fallfac=1.0
integer, parameter :: nadd=5

!-----------------------------------------------------------------------
! version control information
!-----------------------------------------------------------------------

character(len=128)  :: version =  '$Id: zetac_advect_vert.f90,v 1.1.2.4.2.4 2005/06/16 18:01:32 stg Exp $'
character(len=128)  :: tag     =  '$Name:  $'

contains

!#######################################################################

subroutine advect_vert_at_u (                                          &
                                                       Mgrid, Vmetric, &
                                            denu, denm, denh, var, ww, &
                                                                index, &
                                                          advect_vert, &
                                                        do_flux, delt, &
                                                           past, pres, &
                                                               varref )

type (horiz_grid_type),                         intent (in)  ::        &
                                                                Mgrid

type (vert_metric_type),                        intent (in)  ::        &
                                                              Vmetric

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked, ntime),    intent (in)  ::        &
                                                           denu, denm

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),           intent (in)  ::        &
                                                                 denh

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked, ntime),    intent (in)  ::        &
                                                              var, ww

integer,                                        intent (in)  ::        &
                                                                index

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),           intent (out) ::        &
                                                          advect_vert

logical,                                        intent (in)  ::        &
                                                              do_flux

real,                                           intent (in)  ::        &
                                                                 delt

integer,                                        intent (in)  ::        &
                                                           past, pres

real, optional, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,    &
                          Mgrid%kbd:Mgrid%ked), intent (in)  ::        &
                                                               varref

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
     call advect_vert_init ( Mgrid, Vmetric )
  endif

  call get_method ( index, 'advect_vert', name_scheme, num_scheme )
  call get_method ( index, 'equation', name_form, num_form )

  if ( num_scheme == SECOND_CENTERED .and.                             &
       name_scheme(1:3) == "alt" ) then

!-----------------------------------------------------------------------
! second-centered advective form (inline)
!-----------------------------------------------------------------------

     do j=jbc,jec
        do i=ibc,iec
           watuv(i,j,:) = 0.5*(ww(i,j,:,pres)*denm(i,j,:,pres)         &
                           + ww(i+1,j,:,pres)*denm(i+1,j,:,pres))/     &
                                              denu(i,j,:,pres)
        enddo
     enddo

     if ( do_flux ) then

        do k=kbd+1,ked-1
           do j=jbc,jec
              do i=ibc,iec
                 wvar(i,j,k) = watuv(i,j,k)*denh(i,j,k)/dzetaw(k)*     &
                                  (var(i,j,k,pres) + var(i,j,k+1,pres))
              enddo
           enddo
        enddo
        wvar(:,:,kbd) = 0.
        wvar(:,:,ked) = 0.

        do k=kbd+1,ked
           do j=jbc,jec
              do i=ibc,iec
                 advect_vert(i,j,k) = 0.5/denu(i,j,k,pres)*            &
                                          (wvar(i,j,k) - wvar(i,j,k-1))
              enddo
           enddo
        enddo

     else

        do k=kbd+1,ked-1
           do j=jbc,jec
              do i=ibc,iec
                 wvar(i,j,k) = watuv(i,j,k)/dzetaw(k)*                 &
                                  (var(i,j,k+1,pres) - var(i,j,k,pres))
              enddo
           enddo
        enddo
        wvar(:,:,kbd) = 0.
        wvar(:,:,ked) = 0.

        do k=kbd+1,ked
           do j=jbc,jec
              do i=ibc,iec
                 advect_vert(i,j,k) = 0.5*(wvar(i,j,k) + wvar(i,j,k-1))
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
     num_scheme = num_scheme + 5  ! variable grid spacing
     ltime = pres
  else
     ltime = past
  endif

  do k=kbd,ked
     do j=jbc,jec
        do i=ibc,iec
           ww_fms(i,j,k) = 0.5*(ww(i+1,j,k,ltime) + ww(i,j,k,ltime))
        enddo
     enddo
  enddo
  ww_fms(:,:,ked+1) = 0.0

  do k=kbd,ked-1
     do j=jbc,jec
        do i=ibc,iec
           var_fms(i,j,k) = var(i,j,k+1,ltime)
            dz_fms(i,j,k) = dzetam(k+1)
        enddo
     enddo
  enddo
  dz_fms(:,:,ked)  = dz_fms(:,:,ked-1)

  if ( present(varref) ) then
     do j=jbc,jec
        do i=ibc,iec
           var_fms(i,j,ked) = var_fms(i,j,ked-1)                       &
                                + (varref(i,j,ked) - varref(i,j,ked-1))
        enddo
     enddo
  else
     do j=jbc,jec
        do i=ibc,iec
           var_fms(i,j,ked) = var_fms(i,j,ked-1)
        enddo
     enddo
  endif

  if ( do_flux ) then

     do k=kbd,ked
        do j=jbc,jec
           do i=ibc,iec
              dens_fms(i,j,k) = 0.5*(denu(i+1,j,k,ltime)               &
                                   + denu(i  ,j,k,ltime))
           enddo
        enddo
     enddo
     do j=jbc,jec
        do i=ibc,iec
           dens_fms(i,j,ked+1) = dens_fms(i,j,ked)**2/dens_fms(i,j,ked-1)
        enddo
     enddo

     do k=kbd,ked
        do j=jbc,jec
           do i=ibc,iec
              var_fms(i,j,k) = var_fms(i,j,k)*dens_fms(i,j,k+1)
           enddo
        enddo
     enddo

  endif

  dt = delt

  call vert_advection ( dt, ww_fms, dz_fms, var_fms, advect_fms,       &
                        scheme=num_scheme, form=FLUX_FORM )

  do k=kbd+1,ked
     do j=jbc,jec
        do i=ibc,iec
           advect_vert(i,j,k) = -advect_fms(i,j,k-1)
        enddo
     enddo
  enddo

  if ( do_flux ) then

!-----------------------------------------------------------------------
! finish flux-form tendency
!-----------------------------------------------------------------------

     do k=kbd,ked
        do j=jbc,jec
           do i=ibc,iec
              dens_fms(i,j,k) = 0.5*(denu(i+1,j,k,pres)                &
                                   + denu(i  ,j,k,pres))
           enddo
        enddo
     enddo
     do k=kbd+1,ked
        do j=jbc,jec
           do i=ibc,iec
              advect_vert(i,j,k) = advect_vert(i,j,k)/dens_fms(i,j,k)
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
              advect_vert(i,j,k) = advect_vert(i,j,k)                  &
                                                   - var_fms(i,j,k-1)* &
                  (ww_fms(i,j,k) - ww_fms(i,j,k-1)) / dz_fms(i,j,k-1)
           enddo
        enddo
     enddo

  endif

  return
end subroutine advect_vert_at_u

!#######################################################################

subroutine advect_vert_at_v (                                          &
                                                       Mgrid, Vmetric, &
                                            denv, denm, denh, var, ww, &
                                                                index, &
                                                          advect_vert, &
                                                        do_flux, delt, &
                                                           past, pres, &
                                                               varref )

type (horiz_grid_type),                     intent (in)  ::            &
                                                                Mgrid

type (vert_metric_type),                    intent (in)  ::            &
                                                              Vmetric

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked, ntime),    intent (in)  ::        &
                                                           denv, denm

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),           intent (in)  ::        &
                                                                 denh

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked, ntime),    intent (in)  ::        &
                                                              var, ww

integer,                                        intent (in)  ::        &
                                                                index

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),           intent (out) ::        &
                                                          advect_vert

logical,                                        intent (in)  ::        &
                                                              do_flux

real,                                           intent (in)  ::        &
                                                                 delt

integer,                                        intent (in)  ::        &
                                                           past, pres

real, optional, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,    &
                          Mgrid%kbd:Mgrid%ked), intent (in)  ::        &
                                                               varref

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
     call advect_vert_init ( Mgrid, Vmetric )
  endif

  call get_method ( index, 'advect_vert', name_scheme, num_scheme )
  call get_method ( index, 'equation', name_form, num_form )

  if ( num_scheme == SECOND_CENTERED .and.                             &
       name_scheme(1:3) == "alt" ) then

!-----------------------------------------------------------------------
! second-centered advective form (inline)
!-----------------------------------------------------------------------

     do j=jbc,jec
        do i=ibc,iec
           watuv(i,j,:) = 0.5*(ww(i,j,:,pres)*denm(i,j,:,pres)         &
                           + ww(i,j+1,:,pres)*denm(i,j+1,:,pres))/     &
                                              denv(i,j,:,pres)
        enddo
     enddo

     if ( do_flux ) then

        do k=kbd+1,ked-1
           do j=jbc,jec
              do i=ibc,iec
                 wvar(i,j,k) = watuv(i,j,k)*denh(i,j,k)/dzetaw(k)*     &
                                  (var(i,j,k,pres) + var(i,j,k+1,pres))
              enddo
           enddo
        enddo
        wvar(:,:,kbd) = 0.
        wvar(:,:,ked) = 0.

        do k=kbd+1,ked
           do j=jbc,jec
              do i=ibc,iec
                 advect_vert(i,j,k) = 0.5/denv(i,j,k,pres)*            &
                                          (wvar(i,j,k) - wvar(i,j,k-1))
              enddo
           enddo
        enddo

     else

        do k=kbd+1,ked-1
           do j=jbc,jec
              do i=ibc,iec
                 wvar(i,j,k) = watuv(i,j,k)/dzetaw(k)*                 &
                                  (var(i,j,k+1,pres) - var(i,j,k,pres))
              enddo
           enddo
        enddo
        wvar(:,:,kbd) = 0.
        wvar(:,:,ked) = 0.

        do k=kbd+1,ked
           do j=jbc,jec
              do i=ibc,iec
                 advect_vert(i,j,k) = 0.5*(wvar(i,j,k) + wvar(i,j,k-1))
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
     num_scheme = num_scheme + 5  ! variable grid spacing
     ltime = pres
  else 
     ltime = past
  endif

  do k=kbd,ked
     do j=jbc,jec
        do i=ibc,iec
           ww_fms(i,j,k) = (ww(i,j+1,k,ltime) + ww(i,j,k,ltime))*0.5
        enddo
     enddo
  enddo
  ww_fms(:,:,ked+1) = 0.0

  do k=kbd+1,ked
     do j=jbc,jec
        do i=ibc,iec
           var_fms(i,j,k-1) =    var(i,j,k,ltime)
            dz_fms(i,j,k-1) = dzetam(k)
        enddo
     enddo
  enddo
  dz_fms(:,:,ked) =  dz_fms(:,:,ked-1)

  if ( present(varref) ) then
     do j=jbc,jec
        do i=ibc,iec
           var_fms(i,j,ked) = var_fms(i,j,ked-1)                       &
                                + (varref(i,j,ked) - varref(i,j,ked-1))
        enddo
     enddo
  else
     do j=jbc,jec
        do i=ibc,iec
           var_fms(i,j,ked) = var_fms(i,j,ked-1)
        enddo
     enddo
  endif

  if ( do_flux ) then

     do k=kbd,ked
        do j=jbc,jec
           do i=ibc,iec
              dens_fms(i,j,k) = 0.5*(denv(i,j,  k,ltime)               &
                                   + denv(i,j+1,k,ltime))
           enddo
        enddo
     enddo
     do j=jbc,jec
        do i=ibc,iec
           dens_fms(i,j,ked+1) = dens_fms(i,j,ked)**2/dens_fms(i,j,ked-1)
        enddo
     enddo

     do k=kbd,ked
        do j=jbc,jec
           do i=ibc,iec
              var_fms(i,j,k) = var_fms(i,j,k)*dens_fms(i,j,k+1)
           enddo
        enddo
     enddo

  endif

  dt = delt

  call vert_advection ( dt, ww_fms, dz_fms, var_fms, advect_fms,       &
                        scheme=num_scheme, form=FLUX_FORM )

  do k=kbd+1,ked
     do j=jbc,jec
        do i=ibc,iec
           advect_vert(i,j,k) = -advect_fms(i,j,k-1)
        enddo
     enddo
  enddo

  if ( do_flux ) then

!-----------------------------------------------------------------------
! finish flux-form tendency
!-----------------------------------------------------------------------

     do k=kbd,ked
        do j=jbc,jec
           do i=ibc,iec
              dens_fms(i,j,k) = 0.5*(denv(i,j,  k,pres)                &
                                   + denv(i,j+1,k,pres))
           enddo
        enddo
     enddo
     do k=kbd+1,ked
        do j=jbc,jec
           do i=ibc,iec
              advect_vert(i,j,k) = advect_vert(i,j,k)/dens_fms(i,j,k)
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
              advect_vert(i,j,k) = advect_vert(i,j,k)                  &
                                                   - var_fms(i,j,k-1)* &
                  (ww_fms(i,j,k) - ww_fms(i,j,k-1)) / dz_fms(i,j,k-1)
           enddo
        enddo
     enddo

  endif

  return
end subroutine advect_vert_at_v

!#######################################################################

subroutine advect_vert_at_mass (                                       &
                                                       Mgrid, Vmetric, &
                                                  var, ww, denm, denh, &
                                                                index, &
                                                          advect_vert, &
                                                        do_flux, delt, &
                                                           past, pres, &
                                                               varref )

type (horiz_grid_type),                         intent (in)  ::        &
                                                                Mgrid

type (vert_metric_type),                        intent (in)  ::        &
                                                              Vmetric

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked, ntime),    intent (in)  ::        &
                                                        var, ww, denm

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),           intent (in)  ::        &
                                                                 denh

integer,                                        intent (in)  ::        &
                                                                index

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),           intent (out) ::        &
                                                          advect_vert

logical,                                        intent (in)  ::        &
                                                              do_flux

real,                                           intent (in)  ::        &
                                                                 delt

integer,                                        intent (in)  ::        &
                                                           past, pres

real, optional, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,    &
                Mgrid%kbd:Mgrid%ked),           intent (in)  ::        &
                                                               varref

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
     call advect_vert_init ( Mgrid, Vmetric )
  endif

  call get_method ( index, 'advect_vert', name_scheme, num_scheme )
  call get_method ( index, 'equation', name_form, num_form )

  if ( num_scheme == SECOND_CENTERED .and.                             &
       name_scheme(1:3) == "alt" ) then

!-----------------------------------------------------------------------
! second-centered flux or advective form (inline)
!-----------------------------------------------------------------------

     if ( do_flux ) then

        do k=kbd+1,ked-1
           do j=jbc,jec
              do i=ibc,iec
                 wvar(i,j,k) = ww(i,j,k,pres)*denh(i,j,k)/dzetaw(k)*   &
                                 (var(i,j,k,pres) + var(i,j,k+1,pres))
              enddo
           enddo
        enddo
        wvar(:,:,kbd) = 0.
        wvar(:,:,ked) = 0.

        do k=kbd+1,ked
           do j=jbc,jec
              do i=ibc,iec
                  advect_vert(i,j,k) = 0.5/denm(i,j,k,pres)*           &
                                          (wvar(i,j,k) - wvar(i,j,k-1))
              enddo
           enddo
        enddo

     else

        do k=kbd+1,ked-1
           do j=jbc,jec
              do i=ibc,iec
!!$                  wvar(i,j,k) = ww(i,j,k,pres)/dzetaw(k)*              &
!!$                                  (var(i,j,k+1,pres) - var(i,j,k,pres))
                  wvar(i,j,k) = ww(i,j,k,pres)*denh(i,j,k)/dzetaw(k)*  &
                                  (var(i,j,k+1,pres) - var(i,j,k,pres))
              enddo
           enddo
        enddo
        wvar(:,:,kbd) = 0.
        wvar(:,:,ked) = 0.

        do k=kbd+1,ked
           do j=jbc,jec
              do i=ibc,iec
                 advect_vert(i,j,k) = 0.5/denm(i,j,k,pres)*            &
                                          (wvar(i,j,k) + wvar(i,j,k-1))
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
     num_scheme = num_scheme + 5  ! variable grid spacing
     ltime = pres
  else 
     ltime = past
  endif

  do k=kbd,ked
     do j=jbc,jec
        do i=ibc,iec
           ww_fms(i,j,k) = ww(i,j,k,pres)
        enddo
     enddo
  enddo
  ww_fms(:,:,ked+1) = 0.0

  do k=kbd+1,ked
     do j=jbc,jec
        do i=ibc,iec
           var_fms(i,j,k-1) =   var(i,j,k,ltime)
            dz_fms(i,j,k-1) = dzetam(k)
        enddo
     enddo
  enddo
  dz_fms(:,:,ked)  = dz_fms(:,:,ked-1)

  if ( present(varref) ) then
     do j=jbc,jec
        do i=ibc,iec
           var_fms(i,j,ked) = var_fms(i,j,ked-1)                       &
                                + (varref(i,j,ked) - varref(i,j,ked-1))
        enddo
     enddo
  else
     do j=jbc,jec
        do i=ibc,iec
           var_fms(i,j,ked) = var_fms(i,j,ked-1)
        enddo
     enddo
  endif

  if ( do_flux ) then

     do k=kbd,ked
        do j=jbc,jec
           do i=ibc,iec
              dens_fms(i,j,k) = denm(i,j,k,ltime)
           enddo
        enddo
     enddo
     do j=jbc,jec
        do i=ibc,iec
           dens_fms(i,j,ked+1) = denm(i,j,ked,ltime)**2/              &
                                 denm(i,j,ked-1,ltime)
        enddo
     enddo
     do k=kbd,ked
        do j=jbc,jec
           do i=ibc,iec
              var_fms(i,j,k) = var_fms(i,j,k)*dens_fms(i,j,k+1)
           enddo
        enddo
     enddo

  endif

  dt = delt

  call vert_advection ( dt, ww_fms, dz_fms, var_fms, advect_fms,       &
                        scheme=num_scheme, form=FLUX_FORM )

  do k=kbd+1,ked
     do j=jbc,jec
        do i=ibc,iec
           advect_vert(i,j,k) = -advect_fms(i,j,k-1)
        enddo
     enddo
  enddo

  if ( do_flux ) then

!-----------------------------------------------------------------------
! finish flux-form tendency
!-----------------------------------------------------------------------

     do k=kbd+1,ked
        do j=jbc,jec
           do i=ibc,iec
              advect_vert(i,j,k) = advect_vert(i,j,k)/                 &
                                          denm(i,j,k,pres)
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
              advect_vert(i,j,k) = advect_vert(i,j,k)                  &
                                                   - var_fms(i,j,k-1)* &
                  (ww_fms(i,j,k) - ww_fms(i,j,k-1)) / dz_fms(i,j,k-1)
           enddo
        enddo
     enddo

  endif

  return
end subroutine advect_vert_at_mass

!#######################################################################

subroutine advect_vert_at_tracer (                                     &
                                                       Mgrid, Vmetric, &
                                                  var, ww, denm, denh, &
                                                                index, &
                                                          advect_vert, &
                                                                 delt, &
                                                   past, pres, future, &
                                                        fall, fallout )

type (horiz_grid_type),                         intent (in)   ::       &
                                                                Mgrid

type (vert_metric_type),                        intent (in)   ::       &
                                                              Vmetric

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked, ntime),    intent (in)  ::        &
                                                        var, ww, denm

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),           intent (in)  ::        &
                                                                 denh

integer,                                        intent (in)  ::        &
                                                                index

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),           intent (out) ::        &
                                                          advect_vert

real,                                           intent (in)  ::        &
                                                                 delt

integer,                                        intent (in)  ::        &
                                                   past, pres, future

real, optional, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,    &
                          Mgrid%kbd:Mgrid%ked), intent (in)   ::       &
                                                                 fall

real, optional, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed),   &
                                                intent (out) ::        &
                                                              fallout

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

real    :: dt, add
integer :: ltime
integer :: i, j, k

!-----------------------------------------------------------------------
! initialize module
!-----------------------------------------------------------------------

  if (do_init) then
     call advect_vert_init ( Mgrid, Vmetric )
  endif
  
  call get_method ( index, 'advect_vert', name_scheme, num_scheme )
  call get_method ( index, 'equation', name_form, num_form )

!-----------------------------------------------------------------------
! second-centered advective form (inline)
!-----------------------------------------------------------------------

  if ( name_scheme(1:4) == 'none' ) then
     do k=kbd+1,ked
        do j=jbc,jec
           do i=ibc,iec
              advect_vert(i,j,k) = 0.0
           enddo
        enddo
     enddo
     if ( present(fallout) ) then
        do j=jbc,jec
           do i=ibc,iec
              fallout(i,j) = 0.0
           enddo
        enddo
        endif
     return
  endif

  if ( num_scheme == SECOND_CENTERED .and.                             &
       name_scheme(1:3) == "alt" .and. .not. present(fall) ) then

     if ( num_form == ADVECTIVE_FORM ) then

        do k=kbd+1,ked-1
           do j=jbc,jec
              do i=ibc,iec
                  wvar(i,j,k) = ww(i,j,k,pres)/dzetaw(k)*              &
                                  (var(i,j,k+1,pres) - var(i,j,k,pres))
              enddo
           enddo
        enddo
        wvar(:,:,kbd) = 0.
        wvar(:,:,ked) = 0.

        do k=kbd+1,ked
           do j=jbc,jec
              do i=ibc,iec
                 advect_vert(i,j,k) = 0.5*(wvar(i,j,k) + wvar(i,j,k-1))
              enddo
           enddo
        enddo

     else

        do k=kbd+1,ked-1
           do j=jbc,jec
              do i=ibc,iec
                 wvar(i,j,k) = ww(i,j,k,pres)*denh(i,j,k)/dzetaw(k)*   &
                                  (var(i,j,k,pres) + var(i,j,k+1,pres))
              enddo
           enddo
        enddo
        wvar(:,:,kbd) = 0.
        wvar(:,:,ked) = 0.

        do k=kbd+1,ked
           do j=jbc,jec
              do i=ibc,iec
                  advect_vert(i,j,k) = 0.5/denm(i,j,k,pres)*           &
                                          (wvar(i,j,k) - wvar(i,j,k-1))
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
     num_scheme = num_scheme + 5  ! variable grid spacing
     ltime = pres
  else 
     ltime = past
  endif

  do k=kbd,ked
     do j=jbc,jec
        do i=ibc,iec
           ww_fms(i,j,k) = ww(i,j,k,pres)
        enddo
     enddo
   enddo
  ww_fms(:,:,ked+1) = 0.0

  do k=kbd+1,ked
     do j=jbc,jec
        do i=ibc,iec
           var_fms(i,j,k-1) =   var(i,j,k,ltime)
            dz_fms(i,j,k-1) = dzetam(k)
        enddo
     enddo
  enddo
  var_fms(:,:,ked) = var_fms(:,:,ked-1)
  dz_fms(:,:,ked) = dz_fms(:,:,ked-1)  

  if ( present(fall) ) then
     do k=kbd,ked
        do j=jbc,jec
           do i=ibc,iec
              fall_fms(i,j,k) = min( fall(i,j,k),                      &
                                     fallfac*dzetaw(k)/(2.0*delt) )
              ww_fms(i,j,k) = ww_fms(i,j,k) - fall_fms(i,j,k)
           enddo
        enddo
     enddo
  endif

  if ( num_form == FLUX_FORM .or. present(fall) ) then

     do k=kbd,ked
        do j=jbc,jec
           do i=ibc,iec
              dens_fms(i,j,k) = denm(i,j,k,ltime)
           enddo
        enddo
     enddo
     do j=jbc,jec
        do i=ibc,iec
           dens_fms(i,j,ked+1) = denm(i,j,ked,ltime)**2/               &
                                 denm(i,j,ked-1,ltime)
        enddo
     enddo
     do k=kbd,ked
        do j=jbc,jec
           do i=ibc,iec
              var_fms(i,j,k) = var_fms(i,j,k)*dens_fms(i,j,k+1)
           enddo
        enddo
     enddo

  endif

  dt = delt

  call vert_advection ( dt, ww_fms, dz_fms, var_fms, advect_fms,       &
                        scheme=num_scheme, form=FLUX_FORM )

  do k=kbd+1,ked
     do j=jbc,jec
        do i=ibc,iec
           advect_vert(i,j,k) = -advect_fms(i,j,k-1)
        enddo
     enddo
  enddo

!-----------------------------------------------------------------------
! diagnose precipitation (correct for unphysical fallspeed convergence)
!-----------------------------------------------------------------------

  if ( present(fallout) ) then

     do j=jbc,jec
        do i=ibc,iec
           fallout(i,j) = fall_fms(i,j,kbd)*var_fms(i,j,kbd)
        enddo
     enddo

     do k=kbd+1,kbd+nadd
        do j=jbc,jec
           do i=ibc,iec
              add = addfac*var_fms(i,j,k-1)*                           &
                      dim(fall_fms(i,j,k), fall_fms(i,j,k-1))
              fallout(i,j) = fallout(i,j) + add
              advect_vert(i,j,k) = advect_vert(i,j,k) + add/dzetam(k)
           enddo
        enddo
     enddo

  endif

  if ( num_form == FLUX_FORM ) then

!-----------------------------------------------------------------------
! finish flux-form tendency
!-----------------------------------------------------------------------

     do k=kbd+1,ked
        do j=jbc,jec
           do i=ibc,iec
              advect_vert(i,j,k) = advect_vert(i,j,k)/                 &
                                          denm(i,j,k,future)
           enddo
        enddo
     enddo

  else

!-----------------------------------------------------------------------
! change back to advect form
!-----------------------------------------------------------------------

     if ( present(fall) ) then
        do k=kbd,ked
           do j=jbc,jec
              do i=ibc,iec
                 ww_fms(i,j,k) = (ww_fms(i,j,k) + fall_fms(i,j,k))*    &
                              0.5*(dens_fms(i,j,k) + dens_fms(i,j,k+1))
              enddo
           enddo
        enddo
     endif
     do k=kbd+1,ked
        do j=jbc,jec
           do i=ibc,iec
              advect_vert(i,j,k) = advect_vert(i,j,k)                  &
                                                   - var(i,j,k,ltime)* &
                      (ww_fms(i,j,k) - ww_fms(i,j,k-1))/dz_fms(i,j,k-1)
           enddo
        enddo
     enddo
     if ( present(fall) ) then
        do k=kbd+1,ked
           do j=jbc,jec
              do i=ibc,iec
                 advect_vert(i,j,k) = advect_vert(i,j,k)/dens_fms(i,j,k)
              enddo
           enddo
        enddo
     endif

  endif

  return
end subroutine advect_vert_at_tracer

!#######################################################################

subroutine advect_vert_quick (                                         &
                                                       Mgrid, Vmetric, &
                                                  var, ww, denm, denh, &
                                                          advect_vert )

type (horiz_grid_type),                         intent (in)  ::        &
                                                                Mgrid

type (vert_metric_type),                        intent (in)  ::        &
                                                              Vmetric

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),           intent (in)  ::        &
                                                  var, ww, denm, denh

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),           intent (out) ::        &
                                                          advect_vert

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

integer :: i, j, k

!-----------------------------------------------------------------------
! initialize module
!-----------------------------------------------------------------------

  if (do_init) then
     call advect_vert_init ( Mgrid, Vmetric )
  endif

!-----------------------------------------------------------------------
! second-centered advective form
!-----------------------------------------------------------------------

     do k=kbd+1,ked-1
        do j=jbc,jec
           do i=ibc,iec
               wvar(i,j,k) = ww(i,j,k)*denh(i,j,k)/dzetaw(k)*          &
                                       (var(i,j,k+1) - var(i,j,k))
           enddo
        enddo
     enddo
     wvar(:,:,kbd) = 0.
     wvar(:,:,ked) = 0.

     do k=kbd+1,ked
        do j=jbc,jec
           do i=ibc,iec
              advect_vert(i,j,k) = 0.5/denm(i,j,k)*                    &
                                      (wvar(i,j,k) + wvar(i,j,k-1))
           enddo
        enddo
     enddo

 end subroutine advect_vert_quick

!#######################################################################

subroutine advect_vert_init ( Mgrid, Vmetric )

type (horiz_grid_type),                     intent (in)  ::            &
                                                                Mgrid

type (vert_metric_type),                    intent (in)  ::            &
                                                              Vmetric

integer :: ibd, ied
integer :: jbd, jed
                                                           
  ibc = Mgrid%ibc
  iec = Mgrid%iec
  jbc = Mgrid%jbc
  jec = Mgrid%jec

  ibd = Mgrid%ibd
  ied = Mgrid%ied
  jbd = Mgrid%jbd
  jed = Mgrid%jed

  kbd = Mgrid%kbd
  ked = Mgrid%ked

!-----------------------------------------------------------------------
! write version number to logfile
!-----------------------------------------------------------------------

  call write_version_number (version, tag)

!-----------------------------------------------------------------------
! allocate work arrays
!-----------------------------------------------------------------------

  allocate (ww_fms    (ibc:iec,jbc:jec,kbd:ked+1))
  allocate (dens_fms  (ibc:iec,jbc:jec,kbd:ked+1))
  allocate (fall_fms  (ibc:iec,jbc:jec,kbd:ked+1))
  allocate (var_fms   (ibc:iec,jbc:jec,kbd:ked))
  allocate (dz_fms    (ibc:iec,jbc:jec,kbd:ked))
  allocate (advect_fms(ibc:iec,jbc:jec,kbd:ked))

  allocate (watuv     (ibd:ied,jbd:jed,kbd:ked))
  allocate (wvar      (ibd:ied,jbd:jed,kbd:ked))

  allocate (dzetaw(kbd:ked))
  allocate (dzetam(kbd:ked))

  dzetaw  = Vmetric%dzetaw
  dzetam  = Vmetric%dzetam

  do_init = .false.

  return
end subroutine advect_vert_init

!#######################################################################

end module zetac_advect_vert_mod
