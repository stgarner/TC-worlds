module zetac_calc_var_tend_mod

use mpp_mod,            only : mpp_pe, mpp_root_pe, mpp_chksum,        &
                               input_nml_file, mpp_sum, mpp_npes
use fms_mod,            only : file_exist, write_version_number,       &
                               close_file, open_namelist_file,         &
                               check_nml_error, stdlog, stdout
use constants_mod,      only : grav, cp_air, rdgas, hlv, pi
use diag_manager_mod,   only : register_diag_field
use field_manager_mod,  only : model_atmos
use time_manager_mod,   only : time_type
use tracer_manager_mod, only : get_tracer_index, get_tracer_names

use zetac_advect_horiz_mod,      only : advect_horiz_at_mass,          &
                                        advect_horiz_at_u,             &
                                        advect_horiz_at_v
use zetac_advect_vert_mod,       only : advect_vert_at_mass,           &
                                        advect_vert_at_u,              &
                                        advect_vert_at_v
use zetac_axes_mod,              only : axid, gridu, gridv, gridm
use zetac_axis_names_mod,        only : xm_, ym_, zm_,                 &
                                        xu_, yu_, xv_, yv_
use zetac_cgrid_interp_mod,      only : cgrid_interp_mass,             &
                                        cgrid_interp_uv
use zetac_convert_var_mod,       only : get_tv, get_pp, get_phi
use zetac_coriolis_mod,          only : coriolis_at_u,                 &
                                        coriolis_at_v
use zetac_extrap_var_mod,        only : extrap_var_init
use zetac_field_names_mod,       only : prog_names, diag_names,        &
                                        uu_, vv_, oo_, ta_,            &
                                        qv_, gz_, pp_
use zetac_horiz_grid_type_mod,   only : horiz_grid_type
use zetac_horiz_mean_mod,        only : horiz_mean_init
use zetac_horiz_metric_type_mod, only : horiz_metric_type
use zetac_lateral_boundary_mod,  only : lateral_boundary,              &
                                        lateral_signal,                &
                                        lateral_boundary_init
use zetac_mix_horiz_mod,         only : mix_horiz, mix_horiz_init
use zetac_moisture_mod,          only : get_qsat
use zetac_ncdf_io_mod,           only : ncdf_fms_write
use zetac_phys_con_mod,          only : rrat, pref
use zetac_press_grad_mod,        only : press_grad_xy,                 &
                                        get_pgf
use zetac_time_pointers_mod,     only : ntime
use zetac_tracer_mod,            only : get_method
use zetac_update_halos_mod,      only : update_halos
use zetac_vert_metric_type_mod,  only : vert_metric_type
use zetac_vert_diff_mod,         only : vdiff_elements

!-----------------------------------------------------------------------
! public interfaces
!-----------------------------------------------------------------------

implicit none
private
public calc_ps_tend, calc_tq_tend, calc_uv_tend, calc_var_tend_init,   &
       get_dp, get_oo, get_gz, get_divint, divergence

real, allocatable, dimension(:)     :: zetam, zetaw, dzetam, dzetaw
real, allocatable, dimension(:)     :: dxu, dxv, cosulat, cosvlat
real, allocatable, dimension(:,:)   :: preg, uref, vref, tref, qref
real, allocatable, dimension(:,:)   :: uzref, tzref, qzref, tvref
real, allocatable, dimension(:,:)   :: topog
real, allocatable, dimension(:,:)   :: udiv, vdiv
real, allocatable, dimension(:,:,:) :: advect_horiz, advect_vert
real, allocatable, dimension(:,:,:) :: nudge, div
real, allocatable, dimension(:,:,:) :: n_uu, n_vv, n_tq
real, allocatable, dimension(:,:,:) :: pp_half, pp_full
real, allocatable, dimension(:,:,:) :: pgf, coriolis, omega
real, allocatable, dimension(:,:,:) :: divsum, consum
real, allocatable, dimension(:,:,:) :: dpdx, dpdy, dfdx, dfdy
real, allocatable, dimension(:,:,:) :: nudg_u, nudg_v, nudg_m
real, allocatable, dimension(:,:,:) :: gzatu, gzatv, fzatu, fzatv
real, allocatable, dimension(:,:,:) :: dzfull, dzhalf
real, allocatable, dimension(:,:,:) :: alflo, alfhi
real, allocatable, dimension(:,:,:) :: alflou, alflov, alfhiu, alfhiv
real, allocatable, dimension(:,:,:) :: alfupx, alfvpy
real, allocatable, dimension(:,:,:) :: phihu, phihv, phifu, phifv
real, allocatable, dimension(:,:)   :: arefu, arefv
real, allocatable, dimension(:)     :: spng

real, allocatable, dimension(:,:,:) :: udt, vdt

character(len=*),  parameter :: module='zetac_calc_var_tend_mod'
character(len=16), parameter :: model='atmos_mod'
character (len=16) :: name, units
character (len=64) :: longname

real, parameter :: silly=-99999.9, tiny=1.e-20
real, parameter :: gcp=Grav/Cp_air

real    :: dy, ainf
real    :: ptop, pbot, delp
real    :: fsum

integer :: ibc, iec, ibd, ied
integer :: jbc, jec, jbd, jed
integer :: kbc, kec, kbd, ked
integer :: iz

logical :: lhopen

integer :: pindex, uindex, vindex, tindex, qindex
integer :: id_udt_nudge, id_vdt_nudge, id_tdt_nudge, id_qdt_nudge
integer :: id_qdt_nudge_int

! namelist:
logical :: mix_fast=.false.

! namelist:
real :: nudgwidx=0.0, nudgwidy=0.0, nudgwidz=0.0
real :: nudgfacx=0.0, nudgfacy=0.0, nudgfacz=0.0
real :: nudgfactq=0.0, nudgfacuv=0.0
real :: nudgbuff=0.0
real :: spngdepth=0.0, rad_mean=0.e3
real :: spngfacuv=0.0, spngfactq=0.0
real :: divergfac=0.0
real :: diffm_const, diffh_const
real :: cdm_const, cdh_const
real :: z0, zh
real :: Tlapse=1.0, qlapse=0.7
real :: gust=1.
real :: shap=0.
real :: tau_diff=0.

logical :: slip_w=.true., slip_e=.true., slip_s=.true., slip_n=.true.

!-----------------------------------------------------------------------
! version control information
!-----------------------------------------------------------------------

character(len=128)  :: version =  '$Id: zetac_calc_var_tend.f90,v 1.1.2.11.2.14 2005/07/30 02:42:38 stg Exp $'
character(len=128)  :: tag     =  '$Name:  $'

contains

!#######################################################################

subroutine calc_ps_tend (                                              &
                                                                Mgrid, &
                                        ps, uu, vv, oo, dmu, dmv, dpm, &
                                                                qtend, &
                                                                 f_ps, &
                                                past, present, future, &
                                                                 delt )

type (horiz_grid_type),                      intent (in)    ::         &
                                                                Mgrid

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                ntime),                      intent (in)    ::         &
                                                                   ps

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked,ntime),  intent (in)    ::         &
                                                               uu, vv

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked,ntime),  intent (out)   ::         &
                                                                   oo

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),        intent (in)    ::         &
                                                             dmu, dmv

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked,ntime),  intent (in)    ::         &
                                                                  dpm

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),        intent (in)    ::         &
                                                                qtend

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed),             &
                                             intent (out)   ::         &
                                                                 f_ps

real,                                        intent (in)    ::         &
                                                                 delt
     
integer,                                     intent (in)    ::         &
                                                past, present, future

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

integer :: i, j, k

!-----------------------------------------------------------------------
!                          terms in ps equation
!-----------------------------------------------------------------------

  call get_divint (       Mgrid, uu(:,:,:,present), vv(:,:,:,present), &
                                         dmu, dmv, dpm(:,:,:,present), &
                                                   qtend, div, divsum )

  do j=jbc,jec
     do i=ibc,iec
        f_ps(i,j) = -divsum(i,j,ked)
     enddo
  enddo

!-----------------------------------------------------------------------
! diagnose sigma-dot
!-----------------------------------------------------------------------

  call get_oo (Mgrid, zetaw, divsum, ps(:,:,present), oo(:,:,:,present))

  return
end subroutine calc_ps_tend

!#######################################################################

subroutine calc_uv_tend (                                              &
                                        Mgrid, Hmetric, Vmetric, Time, &
                                ps, uu, vv, oo, ta, qv, dpu, dpv, dpm, &
                                                             hpu, hpv, &
                                         dmu, dmv, dmm, hmu, hmv, hmm, &
                                           umean, vmean, uudiv, vvdiv, &
                                                   f_uu, f_vv, gz, fz, &
                             diffu, diffv, rdzmu, rdzwu, rdzmv, rdzwv, &
                                 cdu, cdv, dflxu, dflxv, vsatu, vsatv, &
                                             beta_u, gamma_u, ldiag_u, &
                                             beta_v, gamma_v, ldiag_v, &
               udt_mix, vdt_mix, udt_vdif, vdt_vdif, udt_adv, vdt_adv, &
                                                   udt_hpgf, vdt_hpgf, &
                                       do_flux, past, present, future, &
                                                                 delt )

type (horiz_grid_type),                     intent (in)    ::          &
                                                                Mgrid

type (horiz_metric_type),                   intent (in)    ::          &
                                                              Hmetric

type (vert_metric_type),                    intent (in)    ::          &
                                                              Vmetric

type (time_type),                           intent (in)    ::          &
                                                                 Time

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                ntime),                     intent (in)    ::          &
                                                                   ps

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked,ntime), intent (in)    ::          &
                                    uu, vv, oo, ta, qv, dpu, dpv, dpm

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),       intent (in)    ::          &
                                                             hpu, hpv

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),       intent (in)    ::          &
                                         dmu, dmv, dmm, hmu, hmv, hmm

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),       intent (in)    ::          &
                                           umean, vmean, uudiv, vvdiv

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),       intent (out)   ::          &
                                                   f_uu, f_vv, gz, fz

real, dimension(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed,               &
                Mgrid%kbd:Mgrid%ked),       intent (in)    ::          &
                                                         diffu, diffv

real, dimension(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed,               &
                Mgrid%kbd:Mgrid%ked),       intent (inout) ::          &
                                           rdzmu, rdzwu, rdzmv, rdzwv

real, dimension(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed),              &
                                             intent (in)   ::          &
                                 cdu, cdv, dflxu, dflxv, vsatu, vsatv

real, dimension(Mgrid%ibc:Mgrid%iec, Mgrid%jbc:Mgrid%jec,              &
                Mgrid%kbd+1:Mgrid%ked),      intent (out)  ::          &
                   beta_u, gamma_u, ldiag_u, beta_v, gamma_v, ldiag_v

real, dimension(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed,               &
                Mgrid%kbd:Mgrid%ked),       intent (out)   ::          &
               udt_mix, vdt_mix, udt_vdif, vdt_vdif, udt_adv, vdt_adv

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),       intent (out)   ::          &
                                                   udt_hpgf, vdt_hpgf

real,                                       intent (in)    ::          &
                                                                 delt
 
integer,                                    intent (in)    ::          &
                                                past, present, future

logical, intent (in) :: do_flux

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

integer :: i, j, k

!-----------------------------------------------------------------------
!                          terms in uu equation
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! geopotential and horizontal pressure gradient
!-----------------------------------------------------------------------

  call get_pp ( ps(:,:,present), zetaw, ptop, pp_half )
  call get_pp ( ps(:,:,present), zetam, ptop, pp_full )

  call get_gz ( Mgrid, pp_half, pp_full, ta(:,:,:,present),            &
                                         qv(:,:,:,present), gz, fz )

  call press_grad_xy ( Mgrid, pp_half, dpdx, dpdy, fz, dfdx, dfdy )

!-----------------------------------------------------------------------
! horizontal advection
!-----------------------------------------------------------------------

  call advect_horiz_at_u (                                             &
                                              Mgrid, Hmetric, Vmetric, &
                                                 uu, uu, vv, dpu, dpv, &
                                                               uindex, &
                                                         advect_horiz, &
                                                        do_flux, delt, &
                                                        past, present, &
                                                              udt_adv )

  udt_adv = advect_horiz - udt_adv

!-----------------------------------------------------------------------
! vertical advection
!-----------------------------------------------------------------------

  call advect_vert_at_u (                                              &
                                                       Mgrid, Vmetric, &
                                                dpu, dpm, hpu, uu, oo, &
                                                               uindex, &
                                                          advect_vert, &
                                                        do_flux, delt, &
                                                        past, present )

!-----------------------------------------------------------------------
! coriolis force
!-----------------------------------------------------------------------

  call coriolis_at_u ( Mgrid, dpu(:,:,:,present), dpv(:,:,:,present),  &
                               uu(:,:,:,present), vv(:,:,:,present),   &
                                                            coriolis )

!-----------------------------------------------------------------------
! zonal pressure gradient force
!-----------------------------------------------------------------------

  call get_alpha ( Mgrid, gz, fz, dmu, dmv, dmm )

  call get_pgf ( Mgrid, dpdx, alflou, alfhiu, pgf )

  pgf = pgf + dfdx

  udt_hpgf = -pgf

!-----------------------------------------------------------------------
! horizontal diffusion
!-----------------------------------------------------------------------

  call mix_horiz ( Mgrid, uu(:,:,:,past), udt_mix, uindex,             &
                                                dp=dpu(:,:,:,present) )

!-----------------------------------------------------------------------
! nudging of uu
!-----------------------------------------------------------------------

  call nudg_spng (                                                     &
                                     Mgrid, uref, umean, nudg_u, spng, &
                                                 nudgfacuv, spngfacuv )

  if ( divergfac > 0. ) then
  do k=kbd+1,ked
     do j=jbc,jec
        do i=ibc,iec
           nudge(i,j,k) = nudge(i,j,k) - divergfac*spng(k)*uudiv(i,j,k)
        enddo
     enddo
  enddo
  endif

!-----------------------------------------------------------------------
! total tendency of uu
!-----------------------------------------------------------------------

  do k=kbd+1,ked
     do j=jbc,jec
        do i=ibc,iec
           f_uu(i,j,k) = udt_mix(i,j,k) + nudge(i,j,k)                 &
                       - (advect_horiz(i,j,k) + advect_vert(i,j,k))    &
                                    + coriolis(i,j,k) - pgf(i,j,k)
         enddo
     enddo
  enddo

  udt_vdif = nudge ! initialize for 'step_var'

!-----------------------------------------------------------------------
!                           terms in vv equation
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! horizontal advection
!-----------------------------------------------------------------------

  call advect_horiz_at_v (                                             &
                                              Mgrid, Hmetric, Vmetric, &
                                                 vv, uu, vv, dpu, dpv, &
                                                               vindex, &
                                                         advect_horiz, &
                                                        do_flux, delt, &
                                                        past, present, &
                                                              vdt_adv )

  vdt_adv = advect_horiz - vdt_adv

!-----------------------------------------------------------------------
! vertical advection
!-----------------------------------------------------------------------

  call advect_vert_at_v (                                              &
                                                       Mgrid, Vmetric, &
                                                dpv, dpm, hpv, vv, oo, &
                                                               vindex, &
                                                          advect_vert, &
                                                        do_flux, delt, &
                                                        past, present )

!-----------------------------------------------------------------------
! coriolis force
!-----------------------------------------------------------------------

  call coriolis_at_v ( Mgrid, dpu(:,:,:,present), dpv(:,:,:,present),  &
                               uu(:,:,:,present), vv(:,:,:,present),   &
                                                            coriolis )

!-----------------------------------------------------------------------
! meridional pressure gradient force
!-----------------------------------------------------------------------

  call get_pgf ( Mgrid, dpdy, alflov, alfhiv, pgf )

  pgf = pgf + dfdy

  vdt_hpgf = -pgf

!-----------------------------------------------------------------------
! horizontal diffusion
!-----------------------------------------------------------------------

  call mix_horiz ( Mgrid, vv(:,:,:,past), vdt_mix, vindex,             &
                                                dp=dpv(:,:,:,present) )

!-----------------------------------------------------------------------
! nudging of vv
!-----------------------------------------------------------------------

  call nudg_spng (                                                     &
                                     Mgrid, vref, vmean, nudg_v, spng, &
                                                 nudgfacuv, spngfacuv )

  if ( divergfac > 0. ) then
  do k=kbd+1,ked
     do j=jbc,jec
        do i=ibc,iec
           nudge(i,j,k) = nudge(i,j,k) - divergfac*spng(k)*vvdiv(i,j,k)
        enddo
     enddo
  enddo
  endif

!-----------------------------------------------------------------------
! total tendency of vv
!-----------------------------------------------------------------------

  do k=kbd+1,ked
     do j=jbc,jec
        do i=ibc,iec
           f_vv(i,j,k) = vdt_mix(i,j,k) + nudge(i,j,k)                 &
                       - (advect_horiz(i,j,k) + advect_vert(i,j,k))    &
                                    + coriolis(i,j,k) - pgf(i,j,k)
        enddo
     enddo
  enddo

  vdt_vdif = nudge ! initialize for 'step_var'

!-----------------------------------------------------------------------
! prep vertical diffusion of uu and vv
!-----------------------------------------------------------------------

  call cgrid_interp_uv ( Mgrid, iz, gz, gzatu, gzatv )
  call cgrid_interp_uv ( Mgrid, iz, fz, fzatu, fzatv )

  call vdiff_elements (                                                &
                                                                Mgrid, &
                                                          diffu, n_uu, &
                                                        dpu, dmu, hmu, &
                                             cdu, dflxu, gzatu, fzatu, &
                                         rdzmu, rdzwu, dzfull, dzhalf, &
                                             beta_u, gamma_u, ldiag_u, &
                                             do_flux, present, future, &
                                                                 delt )

  udt_vdif = udt_vdif - f_uu

  do j=jbc,jec
     do i=ibc,iec
        f_uu(i,j,ked) = f_uu(i,j,ked)                                  &
                        - cdu(i,j)*rdzmu(i,j,ked)*uu(i,j,ked,present)* &
                            ( vsatu(i,j) - dflxu(i,j) )*rdzwu(i,j,ked)
     enddo
  enddo

  udt_vdif = udt_vdif + f_uu

  call vdiff_elements (                                                &
                                                                Mgrid, &
                                                          diffv, n_vv, &
                                                        dpv, dmv, hmv, &
                                             cdv, dflxv, gzatv, fzatv, &
                                         rdzmv, rdzwv, dzfull, dzhalf, &
                                             beta_v, gamma_v, ldiag_v, &
                                             do_flux, present, future, &
                                                                 delt )

  vdt_vdif = vdt_vdif - f_vv

  do j=jbc,jec
     do i=ibc,iec
        f_vv(i,j,ked) = f_vv(i,j,ked)                                  &
                        - cdv(i,j)*rdzmv(i,j,ked)*vv(i,j,ked,present)* &
                            ( vsatv(i,j) - dflxv(i,j) )*rdzwv(i,j,ked)
     enddo
  enddo

  vdt_vdif = vdt_vdif + f_vv

!-----------------------------------------------------------------------
! start signal-speed calculation at lateral boundaries
!-----------------------------------------------------------------------

  if ( lhopen ) then

     call lateral_signal (                                             &
                                                                Mgrid, &
                                                 ta(:,:,:,past), tref, &
                                                    nudg_m, nudgfactq, &
                          dpm(:,:,:,present), qv(:,:,:,present), delt )
 
    call lateral_boundary ( Mgrid,                                     &
               uu(:,:,:,present), vv(:,:,:,present), f_uu, f_vv, uref )

  endif

  return
end subroutine calc_uv_tend

!#######################################################################

subroutine calc_tq_tend (                                              &
                                        Mgrid, Hmetric, Vmetric, Time, &
                           uu, vv, oo, ta, qv, sm, dpu, dpv, dpm, hpm, &
                                                         tmean, qmean, &
                                                     om, oh, alfomega, &
                                                           f_ta, f_qv, &
                                                               gz, fz, &
                                 diffh, rdzm, rdzw, cdh, dflxh, vsatm, &
                                                    tsurf, qsurf, dzw, &
                                             beta_h, gamma_h, ldiag_h, &
      tdt_mix, qdt_mix, tdt_vdif, qdt_vdif, tdt_adv, qdt_adv, sdt_adv, &
                                       do_flux, past, present, future, &
                                                                 delt )

type (horiz_grid_type),                         intent (in)    ::      &
                                                                Mgrid

type (horiz_metric_type),                       intent (in)    ::      &
                                                              Hmetric

type (vert_metric_type),                        intent (in)    ::      &
                                                              Vmetric

type (time_type),                               intent (in)    ::      &
                                                                 Time

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked, ntime),    intent (inout) ::      &
                                uu, vv, oo, ta, qv, sm, dpu, dpv, dpm

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),           intent (in)    ::      &
                                                                  hpm

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),           intent (in)    ::      &
                                                         tmean, qmean
            
real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),           intent (out)   ::      &
                                                     om, oh, alfomega

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),           intent (out)   ::      &
                                                           f_ta, f_qv

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),           intent (in)    ::      &
                                                               gz, fz

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),           intent (in)    ::      &
                                                                diffh

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),           intent (inout) ::      &
                                                           rdzm, rdzw

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed),             &
                                                intent (in)    ::      &
                                                    cdh, dflxh, vsatm

real, dimension(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed),              &
                                                intent (in)    ::      &
                                                         tsurf, qsurf

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),           intent (out)   ::      &
                                                                  dzw

real, dimension(Mgrid%ibc:Mgrid%iec, Mgrid%jbc:Mgrid%jec,              &
                Mgrid%kbd+1:Mgrid%ked),         intent (out)   ::      &
                                             beta_h, gamma_h, ldiag_h

real, dimension(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed,               &
                Mgrid%kbd:Mgrid%ked),           intent (out)   ::      &
      tdt_mix, qdt_mix, tdt_vdif, qdt_vdif, tdt_adv, qdt_adv, sdt_adv

real,                                           intent (in)    ::      &
                                                                 delt

integer,                                        intent (in)    ::      &
                                                past, present, future

logical, intent (in) :: do_flux

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

integer :: i, j, k

!-----------------------------------------------------------------------
!                           terms in ta equation
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! horizontal advection
!-----------------------------------------------------------------------

  call advect_horiz_at_mass (                                          &
                                              Mgrid, Hmetric, Vmetric, &
                                            ta, uu, vv, dpu, dpv, dpm, &
                                                               tindex, &
                                                         advect_horiz, &
                                                        do_flux, delt, &
                                                        past, present )

!-----------------------------------------------------------------------
! vertical advection
!-----------------------------------------------------------------------

  call advect_vert_at_mass (                                           &
                                                       Mgrid, Vmetric, &
                                                     ta, oo, dpm, hpm, &
                                                               tindex, &
                                                          advect_vert, &
                                                        do_flux, delt, &
                                                        past, present )

  tdt_adv = -(advect_horiz + advect_vert)

!-----------------------------------------------------------------------
! horizontal diffusion
!-----------------------------------------------------------------------

  call mix_horiz ( Mgrid, ta(:,:,:,past), tdt_mix, tindex,             &
                                                dp=dpm(:,:,:,present) )

!-----------------------------------------------------------------------
! nudging of ta
!-----------------------------------------------------------------------

  call nudg_spng (                                                     &
                                     Mgrid, tref, tmean, nudg_m, spng, &
                                                 nudgfactq, spngfactq )

!-----------------------------------------------------------------------
! alpha*dp/dt
!-----------------------------------------------------------------------

  call simm_burr ( Mgrid,                                              &
           dpu(:,:,:,present), dpv(:,:,:,present), dpm(:,:,:,present), &
            uu(:,:,:,present), vv(:,:,:,present), qv(:,:,:,present),   &
                                         fz, divsum, om, oh, alfomega )

  alfomega = alfomega * (1. + qv(:,:,:,present)) / Cp_air

!-----------------------------------------------------------------------
! total tendency of ta
!-----------------------------------------------------------------------

  do k=kbd+1,ked
     do j=jbc,jec
        do i=ibc,iec
           f_ta(i,j,k) = -(advect_horiz(i,j,k) + advect_vert(i,j,k))   &
                     + alfomega(i,j,k) + tdt_mix(i,j,k) + nudge(i,j,k)
        enddo
     enddo
  enddo

  tdt_vdif = nudge  ! initialize for 'step_var'

!-----------------------------------------------------------------------
! Energy diagnostics:
! sum (u*fu + v*fv - 0.5*(u2 + v2)*fp) = 0, 
!   where fu, fv is only flux convergence and centered diff is used.
! Since sum (fT) = 0 and sum (u*-px + v*-py + alfa*omega) = 0,
!   total RHS is -KE*fp
! If continuous in time, LHS is d(TE)/ dt + KE*fp,
!   hence TE is conserved
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!                           terms in qv equation
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! horizontal advection of qv
!-----------------------------------------------------------------------

  call advect_horiz_at_mass (                                          &
                                              Mgrid, Hmetric, Vmetric, &
                                            qv, uu, vv, dpu, dpv, dpm, &
                                                               qindex, &
                                                         advect_horiz, &
                                                        do_flux, delt, &
                                                        past, present )

!-----------------------------------------------------------------------
! vertical advection of qv
!-----------------------------------------------------------------------

  call advect_vert_at_mass (                                           &
                                                       Mgrid, Vmetric, &
                                                     qv, oo, dpm, hpm, &
                                                               qindex, &
                                                          advect_vert, &
                                                        do_flux, delt, &
                                                        past, present )

  qdt_adv = -(advect_horiz + advect_vert)

!-----------------------------------------------------------------------
! horizontal diffusion of qv
!-----------------------------------------------------------------------

  call mix_horiz ( Mgrid, qv(:,:,:,past), qdt_mix, qindex,             &
                                                dp=dpm(:,:,:,present) )

!-----------------------------------------------------------------------
! nudging of qv
!-----------------------------------------------------------------------

  call nudg_spng (                                                     &
                                     Mgrid, qref, qmean, nudg_m, spng, &
                                                 nudgfactq, spngfactq )

!-----------------------------------------------------------------------
! total tendency of qv
!-----------------------------------------------------------------------

  do k=kbd+1,ked
     do j=jbc,jec
        do i=ibc,iec
           f_qv(i,j,k) = -(advect_horiz(i,j,k) + advect_vert(i,j,k))   &
                                       + qdt_mix(i,j,k) + nudge(i,j,k)
        enddo
     enddo
  enddo

  qdt_vdif = nudge ! initialize for 'step_var'

!-----------------------------------------------------------------------
! moist entropy
!-----------------------------------------------------------------------

sdt_adv=0.
if (.true.) then
  call advect_horiz_at_mass (                                          &
                                              Mgrid, Hmetric, Vmetric, &
                                            sm, uu, vv, dpu, dpv, dpm, &
                                                               qindex, &
                                                         advect_horiz, &
                                                        do_flux, delt, &
                                                        past, present, &
                                                              sdt_adv )

  call advect_vert_at_mass (                                           &
                                                       Mgrid, Vmetric, &
                                                     sm, oo, dpm, hpm, &
                                                               qindex, &
                                                          advect_vert, &
                                                        do_flux, delt, &
                                                        past, present )

  sdt_adv = -(advect_horiz + advect_vert)
endif

!-----------------------------------------------------------------------
! prep vertical diffusion of ta and qv
!-----------------------------------------------------------------------

  call vdiff_elements (                                                &
                                                                Mgrid, &
                                                          diffh, n_tq, &
                                         dpm, dpm(:,:,:,present), hpm, &
                                                   cdh, dflxh, gz, fz, &
                                           rdzm, rdzw, dzfull, dzhalf, &
                                             beta_h, gamma_h, ldiag_h, &
                                             do_flux, present, future, &
                                                                 delt )

  dzw = dzhalf  ! return array to zetac.f90

  tdt_vdif = tdt_vdif - f_ta
  qdt_vdif = qdt_vdif - f_qv

  do k=kbd+1,ked  ! turns T_zz into theta_zz
     do j=jbc,jec
        do i=ibc,iec
           f_ta(i,j,k) = f_ta(i,j,k)                                   &
             + gcp*(diffh(i,j,k-1)*rdzw(i,j,k-1)*dzhalf(i,j,k-1)       &
                  - diffh(i,j,k)*rdzw(i,j,k)*dzhalf(i,j,k))*rdzm(i,j,k)
        enddo
     enddo
  enddo

  do j=jbc,jec
     do i=ibc,iec
        udiv(i,j) = cdh(i,j)*rdzm(i,j,ked)*rdzw(i,j,ked)
        f_ta(i,j,ked) = f_ta(i,j,ked) + udiv(i,j)*                     &
                     (vsatm(i,j)*(tsurf(i,j) - Tlapse*dzhalf(i,j,ked)) &
                       - (vsatm(i,j) - dflxh(i,j))*ta(i,j,ked,present))
        f_qv(i,j,ked) = f_qv(i,j,ked) + udiv(i,j)*                     &
                     (vsatm(i,j)*(qsurf(i,j) - qlapse*dzhalf(i,j,ked)) &
                       - (vsatm(i,j) - dflxh(i,j))*qv(i,j,ked,present))
     enddo
  enddo

  tdt_vdif = tdt_vdif + f_ta
  qdt_vdif = qdt_vdif + f_qv

  return
end subroutine calc_tq_tend

!#######################################################################

subroutine nudg_spng (                                                 &
                                               Mgrid, varref, varmean, & 
                                       nudgefcn, spngefcn, nfac, sfac )

type (horiz_grid_type),                         intent (in)   ::       &
                                                                Mgrid

real, dimension(Mgrid%jbd:Mgrid%jed, Mgrid%iz),                        &
                                                intent (in)   ::       &
                                                               varref

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),                                  &
                                                intent (in)   ::       &
                                                              varmean

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),           intent (in)   ::       &
                                                             nudgefcn

real, dimension(Mgrid%kbd:Mgrid%ked),           intent (in)   ::       &
                                                             spngefcn

real,                                           intent (in)   ::       &
                                                           nfac, sfac

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

integer :: i, k

!-----------------------------------------------------------------------
! nudge towards external field and apply sponge to small scales
!-----------------------------------------------------------------------

  do i=ibd,ied
     nudge(i,:,:) = nudgefcn(i,:,:)*nfac*varref(:,:)
  enddo
  do k=kbd+1,ked
     nudge(:,:,k) = nudge(:,:,k) + spngefcn(k)*sfac*varmean(:,:,k)
  enddo

  return
end subroutine nudg_spng

!#######################################################################

subroutine get_divint (      Mgrid, ucomp, vcomp, dmatu, dmatv, dpatm, &
                                                  qdt, diverg, divint )

!-----------------------------------------------------------------------
! get divergence for pressure tendency
!-----------------------------------------------------------------------

type (horiz_grid_type),                         intent (in)   ::       &
                                                                Mgrid

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),           intent (in)   ::       &
                               ucomp, vcomp, dmatu, dmatv, dpatm, qdt

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),           intent (out)  ::       &
                                                               diverg

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked), optional, intent (out)  ::       &
                                                               divint

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked)                ::         um, vm
integer :: i, j, k

!-----------------------------------------------------------------------
! vertical integral of horizontal divergence
!-----------------------------------------------------------------------

  um = ucomp*dmatu
  vm = vcomp*dmatv

  call divergence ( Mgrid, um, vm, diverg )

  if (.not.present(divint)) return

  divint(:,:,kbd) = 0.
  consum(:,:,kbd) = 0.

  do k=kbd+1,ked
     do j=jbc,jec
        do i=ibc,iec
           divint(i,j,k) = divint(i,j,k-1) + diverg(i,j,k)
           consum(i,j,k) = consum(i,j,k-1) - qdt(i,j,k)*dpatm(i,j,k)
        enddo
     enddo
  enddo

  divint = divint + consum

  return
end subroutine get_divint

!#######################################################################

subroutine divergence ( Mgrid, ucomp, vcomp, diverg )

!-----------------------------------------------------------------------
! vertical integral of linearized mass divergence
!-----------------------------------------------------------------------

type (horiz_grid_type),                         intent (in)    ::      &
                                                                Mgrid

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),           intent (in)    ::      &
                                                         ucomp, vcomp

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),           intent (out)   ::      &
                                                               diverg

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

real    :: dxj, dyj, cosj
integer :: i, j, k

  do k=kbd,ked
     do j=jbc-1,jec
        cosj = cosvlat(j)
        do i=ibc-1,iec
           udiv(i,j) = ucomp(i,j,k)
           vdiv(i,j) = vcomp(i,j,k)*cosj
        enddo
     enddo
     do j=jbc,jec
        dxj = dxu(j)
        dyj = dy*cosulat(j)
        do i=ibc,iec
           diverg(i,j,k) = (udiv(i,j) - udiv(i-1,j))/dxj               &
                         + (vdiv(i,j) - vdiv(i,j-1))/dyj
        enddo
     enddo
  enddo

  return
end subroutine divergence

!#######################################################################

subroutine get_oo ( Mgrid, zeta, divint, ps, oo )

!-----------------------------------------------------------------------
! get divergence for pressure tendency
!-----------------------------------------------------------------------

type (horiz_grid_type),                         intent (in)   ::       &
                                                                Mgrid

real, dimension(Mgrid%iz),                                             &
                                                intent (in)   ::       &
                                                                 zeta

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),           intent (in)   ::       &
                                                               divint

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed),             &
                                                 intent (in)  ::       &
                                                                   ps

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),            intent (out)  ::      &
                                                                   oo

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

integer :: i, j, k

!-----------------------------------------------------------------------
! sigma-dot
!-----------------------------------------------------------------------

  do k=kbd+1,ked-1
     do j=jbc,jec
        do i=ibc,iec
           oo(i,j,k) = (divint(i,j,ked)*zeta(k) - divint(i,j,k))/      &
                                                     (ps(i,j) - ptop)
        enddo
     enddo
  enddo

  oo(:,:,ked) = 0.
  oo(:,:,kbd) = 0.

  return
end subroutine get_oo

!#######################################################################

subroutine simm_burr ( Mgrid, delpu, delpv, delpm, ucomp, vcomp, qvap, &
                                 phim, divint, omega, omeg1, alfomega )

!-----------------------------------------------------------------------
! Simmons-Burridge interpolation of alpha*dp/dt in T equation
!-----------------------------------------------------------------------

type (horiz_grid_type),                      intent (in)    ::         &
                                                                Mgrid
      
real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),        intent (in)    ::         &
                              delpu, delpv, delpm, ucomp, vcomp, qvap, &
                                                         phim, divint
     
real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),        intent (out)   ::         &
                                               omega, omeg1, alfomega
            
!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

integer :: i, j, k

!-----------------------------------------------------------------------
! simmons-burridge interpolation
!-----------------------------------------------------------------------

  do k=kbd+1,ked
     do j=jbc,jec
        do i=ibc-1,iec
           alfupx(i,j,k) = ucomp(i,j,k)*delpu(i,j,k)*                  &
             (dpdx(i,j,k-1)*alfhiu(i,j,k) + dpdx(i,j,k)*alflou(i,j,k))
        enddo
     enddo
     do j=jbc-1,jec
        do i=ibc,iec
           alfvpy(i,j,k) = vcomp(i,j,k)*delpv(i,j,k)*                  &
             (dpdy(i,j,k-1)*alfhiv(i,j,k) + dpdy(i,j,k)*alflov(i,j,k))
        enddo
     enddo
  enddo

  do k=kbd+1,ked
     do j=jbc,jec
        do i=ibc,iec
           alfomega(i,j,k) =                                           &
              0.5*((alfupx(i-1,j,k) + alfupx(i,j,k))                   &
                 + (alfvpy(i,j-1,k) + alfvpy(i,j,k)))/delpm(i,j,k)     &
                  -(alflo(i,j,k)*divint(i,j,k)                         &
                  + alfhi(i,j,k)*divint(i,j,k-1))
        enddo
     enddo
  enddo
  
!-----------------------------------------------------------------------
! pressure velocity at half-levels
!-----------------------------------------------------------------------

  do k=kbd+1,ked-1
     do j=jbc,jec
        do i=ibc-1,iec
           alfupx(i,j,k) = 0.5*(ucomp(i,j,k) + ucomp(i,j,k+1))*dpdx(i,j,k)
        enddo
     enddo
     do j=jbc-1,jec
        do i=ibc,iec
           alfvpy(i,j,k) = 0.5*(vcomp(i,j,k) + vcomp(i,j,k+1))*dpdy(i,j,k)
        enddo
     enddo
  enddo
  do j=jbc,jec
     do i=ibc-1,iec
        alfupx(i,j,kbd) = ucomp(i,j,kbd+1)*dpdx(i,j,kbd)
        alfupx(i,j,ked) = ucomp(i,j,ked)*dpdx(i,j,ked)
     enddo
  enddo
  do j=jbc-1,jec
     do i=ibc,iec
        alfvpy(i,j,kbd) = vcomp(i,j,kbd+1)*dpdy(i,j,kbd)
        alfvpy(i,j,ked) = vcomp(i,j,ked)*dpdy(i,j,ked)
     enddo
  enddo

  do k=kbd,ked
     do j=jbc,jec
        do i=ibc,iec
           omeg1(i,j,k) = 0.5*(alfupx(i-1,j,k) + alfupx(i,j,k)        &
                             + alfvpy(i,j-1,k) + alfvpy(i,j,k))       &
                             - divint(i,j,k)
        enddo
     enddo
  enddo

!-----------------------------------------------------------------------
! pressure velocity at full levels
!-----------------------------------------------------------------------

  do k=kbd+1,ked
     do j=jbc,jec
        do i=ibc-1,iec
           alfupx(i,j,k) = ucomp(i,j,k)*(dpdx(i,j,k-1) + dpdx(i,j,k))
        enddo
     enddo
     do j=jbc-1,jec
        do i=ibc,iec
           alfvpy(i,j,k) = vcomp(i,j,k)*(dpdy(i,j,k-1) + dpdy(i,j,k))
        enddo
     enddo
  enddo

  do k=kbd+1,ked
     do j=jbc,jec
        do i=ibc,iec
           omega(i,j,k) = 0.25*(alfupx(i-1,j,k) + alfupx(i,j,k)        &
                              + alfvpy(i,j-1,k) + alfvpy(i,j,k))       &
                           - 0.5*(divint(i,j,k) + divint(i,j,k-1))
        enddo
     enddo
  enddo
  omega(:,:,kbd) = omega(:,:,kbd+1)

  return
end subroutine simm_burr

!#######################################################################

subroutine get_gz ( Mgrid, phalf, pfull, tabs, qvap, gzhalf, gzfull )

!-----------------------------------------------------------------------
! get geopotential
!-----------------------------------------------------------------------

type (horiz_grid_type),                         intent (in)    ::      &
                                                                Mgrid

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),           intent (in)    ::      &
                                             phalf, pfull, tabs, qvap

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),           intent (out)   ::      &
                                                       gzhalf, gzfull

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

integer :: k

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked)                           ::      &
                                                                 tvrt

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed)      ::      &
                                                                 pbar

!-----------------------------------------------------------------------
! first, half levels ...
!-----------------------------------------------------------------------

  call get_tv ( tabs, qvap, tvrt )

  gzhalf(:,:,ked) = Grav*topog(:,:)

!!$  do k=ked,2,-1
!!$     pbar(:,:) = 0.5*(phalf(:,:,k) + phalf(:,:,k-1))
!!$     gzfull(:,:,k) = gzhalf(:,:,k)                                     &
!!$          + Rdgas*tvrt(:,:,k)/pbar(:,:)*(phalf(:,:,k) - pfull(:,:,k))
!!$     gzhalf(:,:,k-1) = gzfull(:,:,k)                                   &
!!$          + Rdgas*tvrt(:,:,k)/pbar(:,:)*(pfull(:,:,k) - phalf(:,:,k-1))
!!$  enddo

  do k=ked,2,-1
     gzfull(:,:,k) = gzhalf(:,:,k)                                     &
                  + Rdgas*tvrt(:,:,k)*log(phalf(:,:,k)/pfull(:,:,k))
     gzhalf(:,:,k-1) = gzfull(:,:,k)                                   &
                  + Rdgas*tvrt(:,:,k)*log(pfull(:,:,k)/phalf(:,:,k-1))
  enddo 

!!$  do k=ked,2,-1
!!$     pbar(:,:) = phalf(:,:,k) + phalf(:,:,k-1)
!!$     gzfull(:,:,k) = gzhalf(:,:,k)                                     &
!!$          + Rdgas*tvrt(:,:,k)/pbar(:,:)*(phalf(:,:,k) - phalf(:,:,k-1))
!!$     gzhalf(:,:,k-1) = gzfull(:,:,k)                                   &
!!$          + Rdgas*tvrt(:,:,k)/pbar(:,:)*(phalf(:,:,k) - phalf(:,:,k-1))
!!$  enddo

  gzfull(:,:,kbd) =  2.*gzhalf(:,:,kbd) - gzhalf(:,:,kbd+1)

  return
end subroutine get_gz

!#######################################################################

subroutine get_alpha ( Mgrid, phih, phif, delpu, delpv, delpm )

type (horiz_grid_type),                         intent (in)    ::      &
                                                                Mgrid

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),           intent (in)    ::      &
                                      phih, phif, delpu, delpv, delpm

integer :: k

!-----------------------------------------------------------------------
! low and high alpha interpolations
!-----------------------------------------------------------------------

  do k=kbd+1,ked
     alflo(:,:,k) = (phif(:,:,k) - phih(:,:,k))/delpm(:,:,k)
     alfhi(:,:,k) = (phih(:,:,k-1) - phif(:,:,k))/delpm(:,:,k)
  enddo

  call cgrid_interp_uv ( Mgrid, iz, phih, phihu, phihv )
  call cgrid_interp_uv ( Mgrid, iz, phif, phifu, phifv )

  do k=kbd+1,ked
     alflou(:,:,k) = (phifu(:,:,k) - phihu(:,:,k))/delpu(:,:,k)
     alflov(:,:,k) = (phifv(:,:,k) - phihv(:,:,k))/delpv(:,:,k)
     alfhiu(:,:,k) = (phihu(:,:,k-1) - phifu(:,:,k))/delpu(:,:,k)
     alfhiv(:,:,k) = (phihv(:,:,k-1) - phifv(:,:,k))/delpv(:,:,k)
  enddo

  return
end subroutine get_alpha

!#######################################################################

subroutine get_dp ( Mgrid, dzeta, psfc, dpk )

type (horiz_grid_type),                         intent (in)    ::      &
                                                                Mgrid

real, dimension(Mgrid%kbd:Mgrid%ked),                                  &
                                                intent (in)    ::      &
                                                                dzeta

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed),             &
                                                intent (in)    ::      &
                                                                 psfc

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),           intent (out)   ::      &
                                                                  dpk

integer :: k

!-----------------------------------------------------------------------
! pressure thickness dp(k)
!-----------------------------------------------------------------------

  do k=kbd+1,ked
     dpk(:,:,k)= (psfc(:,:) - ptop)*dzeta(k)
  enddo
  dpk(:,:,kbd) = dpk(:,:,kbd+1)

  return
end subroutine get_dp

!#######################################################################

subroutine calc_var_tend_init (                                        &
                                        Mgrid, Hmetric, Vmetric, Time, &
                             ppref, uuref, taref, qvref, gzref, zfref, &
                                                         Tlaps, qlaps, &
                                                  facuv, factq, facsp, &
                                                                 delt, &
                                                             ndt_mean )

type (horiz_grid_type),                     intent (in)    ::          &
                                                                Mgrid

type (horiz_metric_type),                   intent (in)    ::          &
                                                              Hmetric

type (vert_metric_type),                    intent (in)    ::          &
                                                              Vmetric

type (time_type),                           intent (in)    ::          &
                                                                 Time

real, dimension(Mgrid%jbg:Mgrid%jeg, Mgrid%iz),                        &
                                            intent (in)    ::          &
                                           ppref, uuref, taref, qvref

real, dimension(Mgrid%jbg:Mgrid%jeg, Mgrid%iz),                        &
                                            intent (out)   ::          &
                                                                gzref

real, dimension(Mgrid%jbg:Mgrid%jeg),                                  &
                                            intent (out)   ::          &
                                                                zfref

real, intent(out) ::                                     Tlaps, qlaps

real, intent(out) ::                              facuv, factq, facsp

real,                                       intent (in)    ::          &
                                                                  delt

integer,                                    intent (out)    ::         &
                                                             ndt_mean

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

! namelist:
real :: kh2=0.0, kh4=0.0, tqfac=1.0, zonfac=1.0

namelist /mix_horiz_nml/ kh2, kh4, tqfac, zonfac, mix_fast 

namelist /nudging_nml/ nudgwidx, nudgwidy, nudgwidz,                   &
                       nudgfacx, nudgfacy, nudgfacz,                   &
                       nudgfactq, nudgfacuv,                           &
                       nudgbuff,                                       &
                       slip_w, slip_e, slip_s, slip_n

namelist /sponge_nml/ spngdepth, rad_mean,                             &
                      spngfacuv, spngfactq,                            &
                      divergfac,                                       &
                      ndt_mean

namelist /vert_turb_nml/                                               &
                        diffm_const, diffh_const,                      &
                        cdm_const, cdh_const,                          &
                        z0, zh,                                        &
                        Tlapse, qlapse, gust, shap,                    &
                        tau_diff

real :: rlonmin, rlonmax, rlatmin, rlatmax
real :: alon, alat, azee

real, dimension(Mgrid%ibd:Mgrid%ied) :: wtx
real, dimension(Mgrid%jbd:Mgrid%jed) :: wty
real, dimension(Mgrid%iz)            :: wtz
real, dimension(Mgrid%ibd:Mgrid%ied) :: vlon
real, dimension(Mgrid%jbd:Mgrid%jed) :: ulat

real, dimension(Mgrid%jbg:Mgrid%jeg) :: pfref

integer, dimension(2) :: axid_mm
integer, dimension(3) :: axid_uum, axid_vvm, axid_mmm

integer :: io, ierr
integer :: num_form

integer :: i, j, k

  ibc = Mgrid%ibc
  iec = Mgrid%iec
  ibd = Mgrid%ibd
  ied = Mgrid%ied

  jbc = Mgrid%jbc
  jec = Mgrid%jec
  jbd = Mgrid%jbd
  jed = Mgrid%jed

  kbc = Mgrid%kbc
  kec = Mgrid%kec
  kbd = Mgrid%kbd
  ked = Mgrid%ked

  iz = Mgrid%iz
  
  rlonmin = Hmetric%rlonmin
  rlonmax = Hmetric%rlonmax
  rlatmin = Hmetric%rlatmin
  rlatmax = Hmetric%rlatmax
  vlon    = Hmetric%vlon(ibd:ied)
  ulat    = Hmetric%ulat(jbd:jed)

!-----------------------------------------------------------------------
! read namelists
!-----------------------------------------------------------------------

  read (input_nml_file, nml=mix_horiz_nml, iostat=io)
  ierr = check_nml_error(io,'mix_horiz_nml')
  read (input_nml_file, nml=nudging_nml, iostat=io)
  ierr = check_nml_error(io,'nudging_nml')
  read (input_nml_file, nml=sponge_nml, iostat=io)
  ierr = check_nml_error(io,'sponge_nml')
  read (input_nml_file, nml=vert_turb_nml, iostat=io)
  ierr = check_nml_error(io,'vert_turb_nml')

  facuv = spngfacuv
  factq = spngfactq
  facsp = divergfac

!-----------------------------------------------------------------------
! write version number and namelist to logfile
!-----------------------------------------------------------------------

  call write_version_number (version, tag)

  if (mpp_pe() == mpp_root_pe())                                       &
                                    write (stdlog(), nml=mix_horiz_nml)
  if (mpp_pe() == mpp_root_pe())                                       &
                                      write (stdlog(), nml=nudging_nml)
  if (mpp_pe() == mpp_root_pe())                                       &
                                       write (stdlog(), nml=sponge_nml)
 
!-----------------------------------------------------------------------
! allocate and define static arrays
!-----------------------------------------------------------------------

  allocate (zetam(iz), zetaw(iz))
  allocate (dzetam(iz), dzetaw(iz))

  zetam  = Vmetric%zetam
  zetaw  = Vmetric%zetaw
  dzetam  = Vmetric%dzetam
  dzetaw  = Vmetric%dzetaw

  ptop = Vmetric%ptop
  pbot = Vmetric%pbot
  delp = pbot - ptop

  allocate (topog(ibd:ied,jbd:jed))

  topog = Vmetric%topog

  allocate (dxu    (jbd:jed))
  allocate (dxv    (jbd:jed))
  allocate (cosulat(jbd:jed))
  allocate (cosvlat(jbd:jed))

  dy      = Hmetric%dy
  dxu     = Hmetric%dxu(jbd:jed)
  dxv     = Hmetric%dxv(jbd:jed)
  cosulat = Hmetric%cosulat(jbd:jed)
  cosvlat = Hmetric%cosvlat(jbd:jed)

!-----------------------------------------------------------------------
! weights for zonal and meridional sponges
!-----------------------------------------------------------------------

  if ( nudgwidx > 0.0 ) then
     do i=ibd,ied
        alon = min(vlon(i) - rlonmin, rlonmax - vlon(i))
        wtx(i) = nudg_func(alon/nudgwidx)
     enddo
  else
     wtx =0.0
  endif
  
  if ( nudgwidy > 0.0 ) then
     do j=jbd,jed
        alat = min(ulat(j) - rlatmin, rlatmax - ulat(j))
        wty(j) = nudg_func(alat/nudgwidy)
     enddo
  else
     wty = 0.0
  endif

  if ( nudgwidz > 0.0 ) then
     do k=kbd+1,ked
        azee = zetam(k)*delp
        wtz(k) = nudg_func(azee/nudgwidz)
     enddo
  else
     wtz = 0.0
  endif

  allocate (nudg_u(ibd:ied,jbd:jed,kbd:ked))
  allocate (nudg_v(ibd:ied,jbd:jed,kbd:ked))
  allocate (nudg_m(ibd:ied,jbd:jed,kbd:ked))
  nudg_m =0.

  do k=kbd+1,ked
if (zetam(k) > .25) exit
     do j=jbd,jed
        do i=ibd,ied
           nudg_m(i,j,k) = max( wtx(i)*nudgfacx, wty(j)*nudgfacy,      &
                                                 wtz(k)*nudgfacz )
        enddo
     enddo
  enddo

  call cgrid_interp_uv ( Mgrid, iz, nudg_m, nudg_u, nudg_v )
  call update_halos ( Mgrid, nudg_u )
  call update_halos ( Mgrid, nudg_v )

  allocate (spng(kbd:ked))
  spng = 0.

  if ( spngdepth > 0.0 ) then
     do k=kbd+1,ked
if (zetam(k) > .25) exit
        azee = zetam(k)*delp
        spng(k) = spng_func(azee/spngdepth)
     enddo
  endif

if (mpp_pe() == 0) print *,"nudge", nudg_m(Mgrid%ibc,Mgrid%jbc,kbd+1:ked)
if (mpp_pe() == 0) print *,"sponge", spng(kbd+1:ked)
  allocate (n_uu(ibd:ied, jbd:jed, kbd:ked))
  allocate (n_vv(ibd:ied, jbd:jed, kbd:ked))
  allocate (n_tq(ibd:ied, jbd:jed, kbd:ked))

  do k=kbd+1,ked
     do j=jbc,jec
        do i=ibc,iec
           n_uu(i,j,k) = nudgfacuv*nudg_u(i,j,k) + spngfacuv*spng(k)
           n_vv(i,j,k) = nudgfacuv*nudg_v(i,j,k) + spngfacuv*spng(k)
           n_tq(i,j,k) = nudgfactq*nudg_m(i,j,k) + spngfactq*spng(k)
       enddo
     enddo
  enddo

!-----------------------------------------------------------------------
! reference fields
!-----------------------------------------------------------------------

  Tlapse = Tlapse*gcp
  qlapse = qlapse*1.e-5

  Tlaps = Tlapse  ! out to zetac.f90 for diagnostics
  qlaps = qlapse

  allocate (tvref(Mgrid%jbg:Mgrid%jeg, Mgrid%iz))

  call get_tv ( taref, qvref, tvref )
  call get_phi ( ppref, tvref, gzref )
  pfref(:) = ptop + (ppref(:,iz) - ptop)*zetam(iz)
  zfref(:) = Rdgas*tvref(:,iz)/Grav*log(ppref(:,iz)/pfref(:))

  allocate (preg(jbd:jed,kbd:ked))
  allocate (uref(jbd:jed,kbd:ked))
  allocate (vref(jbd:jed,kbd:ked))
  allocate (tref(jbd:jed,kbd:ked))
  allocate (qref(jbd:jed,kbd:ked))

  preg = ppref(jbd:jed,:)
  uref = uuref(jbd:jed,:)
  tref = taref(jbd:jed,:)
  qref = qvref(jbd:jed,:)
  vref = 0.0

  allocate (uzref(jbd:jed,kbd:ked))
  allocate (tzref(jbd:jed,kbd:ked))
  allocate (qzref(jbd:jed,kbd:ked))

  do k=kbd+1,ked-1
     do j=jbd,jed
        uzref(j,k) = (uref(j,k+1) - uref(j,k))/dzetaw(k)
        tzref(j,k) = (tref(j,k+1) - tref(j,k))/dzetaw(k)
        qzref(j,k) = (qref(j,k+1) - qref(j,k))/dzetaw(k)
     enddo
  enddo
  do j=jbd,jed
     uzref(j,kbd) = 0.
     uzref(j,ked) = 0.
     tzref(j,kbd) = 0.
     tzref(j,ked) = 0
     qzref(j,kbd) = 0.
     qzref(j,ked) = 0.
  enddo

  allocate (arefu(jbd:jed, kbd:ked))
  allocate (arefv(jbd:jed, kbd:ked))

  do k=kbd+1,ked
     do j=jbc,jec+1
        arefu(j,k) = Rdgas*tref(j,k)*log(preg(j,k)/preg(j,k-1))/       &
                                        (preg(j,k) - preg(j,k-1))
      enddo
  enddo

  do j=jbc,jec
     arefv(j,:) = 0.5*(arefu(j,:) + arefu(j+1,:))
  enddo

!-----------------------------------------------------------------------
! allocate work arrays
!-----------------------------------------------------------------------

  allocate (advect_horiz(ibd:ied, jbd:jed, kbd:ked))
  allocate (advect_vert (ibd:ied, jbd:jed, kbd:ked))
  allocate (coriolis    (ibd:ied, jbd:jed, kbd:ked))
  allocate (omega       (ibd:ied, jbd:jed, kbd:ked))
  allocate (nudge       (ibd:ied, jbd:jed, kbd:ked))
  allocate (div         (ibd:ied, jbd:jed, kbd:ked))
  allocate (pgf         (ibd:ied, jbd:jed, kbd:ked))
  allocate (dpdx        (ibd:ied, jbd:jed, kbd:ked))
  allocate (dpdy        (ibd:ied, jbd:jed, kbd:ked))
  allocate (dfdx        (ibd:ied, jbd:jed, kbd:ked))
  allocate (dfdy        (ibd:ied, jbd:jed, kbd:ked))
  allocate (divsum      (ibd:ied, jbd:jed, kbd:ked))
  allocate (consum      (ibd:ied, jbd:jed, kbd:ked))
  allocate (alflo       (ibd:ied, jbd:jed, kbd:ked))
  allocate (alfhi       (ibd:ied, jbd:jed, kbd:ked))
  allocate (alflou      (ibd:ied, jbd:jed, kbd:ked))
  allocate (alfhiu      (ibd:ied, jbd:jed, kbd:ked))
  allocate (alflov      (ibd:ied, jbd:jed, kbd:ked))
  allocate (alfhiv      (ibd:ied, jbd:jed, kbd:ked))
  allocate (alfupx      (ibd:ied, jbd:jed, kbd:ked))
  allocate (alfvpy      (ibd:ied, jbd:jed, kbd:ked))
  allocate (pp_half     (ibd:ied, jbd:jed, kbd:ked))
  allocate (pp_full     (ibd:ied, jbd:jed, kbd:ked))
  allocate (phihu       (ibd:ied, jbd:jed, kbd:ked))
  allocate (phifu       (ibd:ied, jbd:jed, kbd:ked))
  allocate (phihv       (ibd:ied, jbd:jed, kbd:ked))
  allocate (phifv       (ibd:ied, jbd:jed, kbd:ked))

  allocate (udiv(ibd:ied, jbd:jed))
  allocate (vdiv(ibd:ied, jbd:jed))

  allocate (gzatu(ibd:ied, jbd:jed, kbd:ked))
  allocate (gzatv(ibd:ied, jbd:jed, kbd:ked))
  allocate (fzatu(ibd:ied, jbd:jed, kbd:ked))
  allocate (fzatv(ibd:ied, jbd:jed, kbd:ked))

  allocate (dzfull(ibd:ied, jbd:jed, kbd:ked))
  allocate (dzhalf(ibd:ied, jbd:jed, kbd:ked))

  allocate (udt(ibd:ied, jbd:jed, kbd:ked))
  allocate (vdt(ibd:ied, jbd:jed, kbd:ked))

  advect_horiz = 0.0
  advect_vert  = 0.0
  coriolis = 0.0
  omega = 0.0
  nudge = 0.0
  pgf = 0.0
  dpdx = 0.0 ; dpdy = 0.0
  dfdx = 0.0 ; dfdy = 0.0
  divsum = 0.0
  alflo = 1.0
  alfhi = 1.0
  alflou = 1.0
  alflov = 1.0
  alfhiu = 1.0
  alfhiv = 1.0

!-----------------------------------------------------------------------
! initialize horizontal diffusion
!-----------------------------------------------------------------------

  call mix_horiz_init (                                                &
                                                       Mgrid, Hmetric, &
                                                             kh2, kh4, &
                                                        tqfac, zonfac, &
                                                                 delt )

!-----------------------------------------------------------------------
! initialize horizontal averaging
!-----------------------------------------------------------------------

  call horiz_mean_init ( Mgrid, Hmetric, rad_mean, facsp, spng )

!-----------------------------------------------------------------------
! prepare extrapolations at open lateral boundaries
!-----------------------------------------------------------------------

  lhopen = ( Mgrid%lxopen .or. Mgrid%lyopen )

  call extrap_var_init (                                               &
                                             Mgrid, Vmetric, nudgbuff, &
                                       slip_w, slip_e, slip_s, slip_n, &
                                                               lhopen )

  if ( lhopen ) then
     call lateral_boundary_init ( Mgrid, Hmetric, Vmetric, delt )
  endif

  pindex = get_tracer_index ( model_atmos, diag_names(pp_) )
  uindex = get_tracer_index ( model_atmos, prog_names(uu_) )
  vindex = get_tracer_index ( model_atmos, prog_names(vv_) )
  tindex = get_tracer_index ( model_atmos, prog_names(ta_) )
  qindex = get_tracer_index ( model_atmos, prog_names(qv_) )

!-----------------------------------------------------------------------
! register diagnostics
!-----------------------------------------------------------------------

  axid_uum = (/axid(xu_),axid(yu_),axid(zm_)/)
  axid_vvm = (/axid(xv_),axid(yv_),axid(zm_)/)
  axid_mmm = (/axid(xm_),axid(ym_),axid(zm_)/)
  axid_mm  = (/axid(xm_),axid(ym_)/)

!!$  call get_tracer_names ( model_atmos, uindex, name, longname, units )
!!$  name = name(1:1) // 'dt_nudge'
!!$  longname = 'nudging tendency of ' // trim(longname)
!!$  units = units // '/s'
!!$  id_udt_nudge = register_diag_field (                                 &
!!$                      model, name, axid_uum, Time,                     &
!!$                      longname, units, silly, range=(/-1.e-2, 1.e-2/) )
!!$
!!$  call get_tracer_names ( model_atmos, vindex, name, longname, units )
!!$  name = name(1:1) // 'dt_nudge'
!!$  longname = 'nudging tendency of ' // trim(longname)
!!$  units = units // '/s'
!!$  id_vdt_nudge = register_diag_field (                                 &
!!$                      model, name, axid_vvm, Time,                     &
!!$                      longname, units, silly, range=(/-1.e-2, 1.e-2/) )
!!$
!!$  call get_tracer_names ( model_atmos, tindex, name, longname, units )
!!$  name = name(1:1) // 'dt_nudge'
!!$  longname = 'nudging tendency of ' // trim(longname)
!!$  units = units // '/s'
!!$  id_tdt_nudge = register_diag_field (                                 &
!!$                       model, name, axid_mmm, Time,                    &
!!$                       longname, units, silly, range=(/-1.e-2,1.e-2/) )
!!$
  call get_tracer_names ( model_atmos, qindex, name, longname, units )
  name = name(1:1) // 'dt_nudge'
  longname = 'nudging tendency of ' // trim(longname)
  units = units // '/s'
  id_qdt_nudge = register_diag_field (                                 &
                       model, name, axid_mmm, Time,                    &
                       longname, units, silly, range=(/-1.e-4,1.e-4/) )
!!$
!!$  name = trim(name) // '_int'
!!$  longname = 'vertically integrated ' // trim(longname)
!!$  units = 'kg/m2/s'
!!$  id_qdt_nudge_int = register_diag_field (                             &
!!$                       model, name, axid_mm, Time,                     &
!!$                       longname, units, silly, range=(/-1.e-2,1.e-2/) )

  ainf = Mgrid%ix*sum(Hmetric%cosulat(1:Mgrid%iy))

  return
end subroutine calc_var_tend_init

!#######################################################################

function nudg_func ( arg )

real, intent (in)  :: arg
real               :: nudg_func

!  nudg_func = 0.5 + 0.5*cos(pi*min(1.0, arg))
!  nudg_func = .5 + sign(.5,1.-arg)
!  nudg_func = dim(1.0, arg)
  nudg_func = 1.

end function nudg_func

!#######################################################################

function spng_func ( arg )

real, intent (in)  :: arg
real               :: spng_func

!  spng_func = 0.5 + 0.5*cos(pi*min(1.0, arg))
!  spng_func = .5 + sign(.5,1.-arg)
!  spng_func = dim(1.0, arg)
  spng_func = 1.

end function spng_func

!#######################################################################

end module zetac_calc_var_tend_mod
