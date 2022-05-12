module zetac_mod

!-----------------------------------------------------------------------
! module information
!-----------------------------------------------------------------------

  use mpp_mod,                   only : mpp_pe, mpp_root_pe,           &
                                        mpp_sum, mpp_npes,             &
                                        mpp_clock_id,                  &
					mpp_clock_begin, mpp_clock_end
use mpp_io_mod,                  only : axistype, fieldtype,           &
                                        mpp_write, mpp_write_meta
use fms_mod,                     only : write_version_number, stdout
use constants_mod,               only : Grav,Cp_air, Hlv
use diag_manager_mod,            only : register_diag_field,           &
                                        need_data
use field_manager_mod,           only : model_atmos
use time_manager_mod,            only : time_type, get_time,           &
                                        operator(-), operator(.eq.)
use tracer_manager_mod,          only : get_tracer_index
use vert_advection_mod,          only : FLUX_FORM

use zetac_axes_mod,              only : axid, gridm
use zetac_axis_names_mod,        only : xm_, ym_, zm_
use zetac_calc_var_tend_mod,     only : calc_ps_tend,                  &
                                        calc_uv_tend,                  &
                                        calc_tq_tend,                  &
                                        calc_var_tend_init,            &
                                        get_dp, get_oo, get_divint,    &
                                        divergence
use zetac_cgrid_interp_mod,      only : cgrid_interp_uv,               &
                                        cgrid_interp_mass
use zetac_convert_var_mod,       only : get_pp, get_ro, get_rh, get_sm
use zetac_coriolis_mod,          only : coriolis_init
use zetac_field_names_mod,       only : prog_names,                    &
                                        uu_, vv_, ta_, qv_, qc_
use zetac_horiz_grid_type_mod,   only : horiz_grid_type
use zetac_horiz_mean_mod,        only : horiz_mean
use zetac_horiz_metric_type_mod, only : horiz_metric_type
use zetac_moisture_mod,          only : get_qsat
use zetac_ncdf_io_mod,           only : ncdf_fms_write,                &
                                        ncdf_fms_write_integral
use zetac_phys_con_mod,          only : rrat
use zetac_press_grad_mod,        only : press_grad_init
use zetac_step_var_mod,          only : step_ps, step_tq, step_uv,     &
                                        step_var_init
use zetac_time_filter_mod,       only : time_filter_init
use zetac_time_pointers_mod,     only : get_time_pointers, ntime
use zetac_tracer_mod,            only : get_method
use zetac_update_halos_mod,      only : update_halos
use zetac_vert_metric_type_mod,  only : vert_metric_type
use zetac_vert_diff_mod,         only : vert_diff_init, vdiff_eval,    &
                                        vdifk_eval
use zetac_vert_turb_mod,         only : vert_turb_init

!-----------------------------------------------------------------------
! public interfaces
!-----------------------------------------------------------------------

implicit none
private
public zetac, zetac_init, zetac_end

character(len=*),  parameter :: module='zetac_mod'
character(len=16), parameter :: model='atmos_mod'
character (len=16) :: name, units
character (len=64) :: longname

integer :: id_udt_vdif, id_vdt_vdif
integer :: id_udt_damp, id_vdt_damp
integer :: id_udt_vdif_int, id_vdt_vdif_int
integer :: id_udt_damp_int, id_vdt_damp_int
integer :: id_divdt_hpgf, id_pgfx, id_pgfy
integer :: id_tsurf, id_qsurf, id_hsurf

real, parameter :: silly=-99999.9
real, parameter :: gcp = Grav/Cp_air
real :: Tlapse, qlapse

real, allocatable, dimension(:) :: zetam, zetaw, dzetam, dzetaw

!-----------------------------------------------------------------------
! sea-surface and air-surface fields
!-----------------------------------------------------------------------

real, allocatable, dimension(:,:) :: vbdy, dzbot
real, allocatable, dimension(:,:) :: tsurf, qsurf, hsurf

!-----------------------------------------------------------------------
! explicit tendencies
!-----------------------------------------------------------------------

real, allocatable, dimension(:,:)   :: f_ps
real, allocatable, dimension(:,:,:) :: f_uu
real, allocatable, dimension(:,:,:) :: f_vv
real, allocatable, dimension(:,:,:) :: f_ta
real, allocatable, dimension(:,:,:) :: f_qv

!-----------------------------------------------------------------------
! horizontal pressure gradient
!-----------------------------------------------------------------------

real, allocatable, dimension(:,:,:) :: udt_hpgf, vdt_hpgf

!-----------------------------------------------------------------------
! density and virtual potential temperature
!-----------------------------------------------------------------------

real, allocatable, dimension(:,:,:) :: beta_u
real, allocatable, dimension(:,:,:) :: beta_v
real, allocatable, dimension(:,:,:) :: beta_h
real, allocatable, dimension(:,:,:) :: gamma_u
real, allocatable, dimension(:,:,:) :: gamma_v
real, allocatable, dimension(:,:,:) :: gamma_h
real, allocatable, dimension(:,:,:) :: ldiag_u
real, allocatable, dimension(:,:,:) :: ldiag_v
real, allocatable, dimension(:,:,:) :: ldiag_h

!-----------------------------------------------------------------------
! interpolation
!----------------------------------------- ------------------------------

real, allocatable, dimension(:,:,:) :: uatm, vatm
real, allocatable, dimension(:,:,:) :: utmp, vtmp
real, allocatable, dimension(:,:,:) :: kpdt

!-----------------------------------------------------------------------
! density, geopotential and pressure
!----------------------------------------- ------------------------------

real, allocatable, dimension(:,:,:) :: rdzm, rdzw
real, allocatable, dimension(:,:,:) :: rdzmu, rdzmv, rdzwu, rdzwv
real, allocatable, dimension(:,:,:) :: dzw, qvh

type (time_type),save :: Time_base

integer :: zetac_core_driver_clock
integer :: zetac_slow_clock
integer :: zetac_step_clock
integer :: ibc, iec, ibd, ied
integer :: jbc, jec, jbd, jed
integer :: j, k, kbd, ked, iz

!namelists:

logical :: uu_is_fluxed=.false., vv_is_fluxed=.false.,                 &
           ta_is_fluxed=.false., qv_is_fluxed=.false.,                 &
           qc_is_fluxed=.false.,                                       &
           uv_is_fluxed = .false., tq_is_fluxed = .false.
integer :: past, pres, future
integer :: secs
integer :: nsteps_leap, nsteps_forward, ndt_mean, dt_mean
real    :: ptop
real    :: delt
real    :: facuv, factq, facsp

logical :: do_first=.true.

!-----------------------------------------------------------------------
! version control information
!-----------------------------------------------------------------------

character(len=128)  :: version =  '$Id: zetac.f90,v 1.1.2.10.2.10 2005/07/25 21:08:11 stg Ppp $'
character(len=128)  :: tag     =  '$Name:  $'

contains

!#######################################################################

subroutine zetac (                                                     &
                                              Mgrid, Hmetric, Vmetric, &
                                           Time_prev, Time, Time_next, &
                        ps, uu, vv, oo, ta, qv, qc, sm, dpu, dpv, dpm, &
                                                        hpu, hpv, hpm, &
                                         dmu, dmv, dmm, hmu, hmv, hmm, &
                         om, oh, alfomega, um, vm, rh, gz, fz, th, qh, &
                                             ttend, qtend, pedt, kedt, &
                               udt_vdif, vdt_vdif, tdt_vdif, qdt_vdif, &
                               udt_damp, vdt_damp, tdt_damp, qdt_damp, &
                                   udt_mix, vdt_mix, tdt_mix, qdt_mix, &
                               udt_filt, vdt_filt, tdt_filt, qdt_filt, &
                          udt_adv, vdt_adv, tdt_adv, qdt_adv, sdt_adv, &
                               umean, vmean, tmean, qmean, udiv, vdiv, &
                                                  diffu, diffv, diffh, &
    cdu, cdv, cdh, dflxu, dflxv, dflxh, vsatu, vsatv, vsatm, sst, ssq )

!-----------------------------------------------------------------------
! define arguments
!-----------------------------------------------------------------------

type (horiz_grid_type)                                         ::      &
                                                                Mgrid

type (horiz_metric_type),                       intent (in)    ::      &
                                                              Hmetric

type (vert_metric_type),                        intent (in)    ::      &
                                                              Vmetric

type (time_type),                               intent (in)    ::      &
                                           Time_prev, Time, Time_next

real, dimension(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed,               &
                ntime),                         intent (inout) ::      &
                                                                   ps

real, dimension(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed,               &
                Mgrid%kbd:Mgrid%ked, ntime),    intent (inout) ::      &
                            uu, vv, oo, ta, qv, qc, sm, dpu, dpv, dpm

real, dimension(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed,               &
                Mgrid%kbd:Mgrid%ked),           intent (inout) ::      &
                                                        hpu, hpv, hpm

real, dimension(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed,               &
                Mgrid%kbd:Mgrid%ked),           intent (inout) ::      &
                                         dmu, dmv, dmm, hmu, hmv, hmm

real, dimension(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed,               &
                Mgrid%kbd:Mgrid%ked),           intent (out)   ::      &
                         om, oh, alfomega, um, vm, rh, gz, fz, th, qh

real, dimension(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed,               &
                Mgrid%kbd:Mgrid%ked),           intent (inout) ::      &
                                             ttend, qtend, pedt, kedt

real, dimension(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed,               &
                Mgrid%kbd:Mgrid%ked),           intent (out)   ::      &
                               udt_vdif, vdt_vdif, tdt_vdif, qdt_vdif

real, dimension(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed,               &
                Mgrid%kbd:Mgrid%ked),           intent (out)   ::      &
                               udt_damp, vdt_damp, tdt_damp, qdt_damp

real, dimension(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed,               &
                Mgrid%kbd:Mgrid%ked),           intent (inout) ::      &
                                   udt_mix, vdt_mix, tdt_mix, qdt_mix

real, dimension(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed,               &
                Mgrid%kbd:Mgrid%ked),           intent (inout) ::      &
                               udt_filt, vdt_filt, tdt_filt, qdt_filt

real, dimension(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed,               &
                Mgrid%kbd:Mgrid%ked),           intent (out)   ::      &
                          udt_adv, vdt_adv, tdt_adv, qdt_adv, sdt_adv

real, dimension(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed,               &
                Mgrid%kbd:Mgrid%ked),           intent (inout) ::      &
                               umean, vmean, tmean, qmean, udiv, vdiv

real, dimension(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed,               &
                Mgrid%kbd:Mgrid%ked),           intent (out)  ::       &
                                                  diffu, diffv, diffh

real, dimension(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed),              &
                                                intent (in)   ::       &
              cdu, cdv, cdh, dflxu, dflxv, dflxh, vsatu, vsatv, vsatm

real, dimension(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed),              &
                                                intent (in)   ::       &
                                                                  sst

real, dimension(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed),              &
                                                intent (out)  ::       &
                                                                  ssq

!-----------------------------------------------------------------------
! get current time, timesteps, and time indices
!-----------------------------------------------------------------------

  call mpp_clock_begin ( zetac_core_driver_clock )

  call get_time ( Time_next - Time_prev, secs )
  delt = real(secs)

  call get_time_pointers ( past, pres, future )

!-----------------------------------------------------------------------
! first time step
!-----------------------------------------------------------------------

  if (do_first) then
     call cgrid_interp_uv ( Mgrid, iz, dpm(:,:,:,past),                &
                      dpu(:,:,:,past), dpv(:,:,:,past) )
     call cgrid_interp_uv ( Mgrid, iz, dpm(:,:,:,pres),                &
                      dpu(:,:,:,pres), dpv(:,:,:,pres) )
     call update_halos ( Mgrid, dpu(:,:,:,past) )
     call update_halos ( Mgrid, dpv(:,:,:,past) )
     call update_halos ( Mgrid, dpu(:,:,:,pres) )
     call update_halos ( Mgrid, dpv(:,:,:,pres) )
     do_first = .false.
     call get_pp ( ps(:,:,past), zetam, ptop, uatm )
     call get_rh ( uatm, ta(:,:,:,past), qv(:,:,:,past), rh, 1.e-6 )
     call get_sm (ta(:,:,:,past), uatm, qv(:,:,:,past), rh, sm(:,:,:,past))
     call get_pp ( ps(:,:,pres), zetam, ptop, vatm )
     call get_rh ( vatm, ta(:,:,:,pres), qv(:,:,:,pres), rh, 1.e-6 )
     call get_sm (ta(:,:,:,pres), vatm, qv(:,:,:,pres), rh, sm(:,:,:,pres))
  endif

!-----------------------------------------------------------------------
! sea-surface specific humidity
!-----------------------------------------------------------------------

  call get_qsat ( ps(:,:,pres), sst, ssq )

!-----------------------------------------------------------------------
! mass at half-levels
!-----------------------------------------------------------------------

  qh(:,:,kbd+1:ked-1) = 0.5*(qv(:,:,kbd+1:ked-1,pres)                 &
                           + qv(:,:,kbd+2:ked,pres))
  th(:,:,kbd+1:ked-1) = 0.5*(ta(:,:,kbd+1:ked-1,pres)                 &
                           + ta(:,:,kbd+2:ked,pres))
  qh(:,:,kbd) = qv(:,:,kbd,pres)
  qh(:,:,ked) = qv(:,:,ked,pres)
  th(:,:,kbd) = ta(:,:,kbd,pres)
  th(:,:,ked) = ta(:,:,ked,pres)

  call dpdry ( Mgrid, dzetaw, ps(:,:,pres), qh, hpu, hpv, hpm )

  call mpp_clock_begin ( zetac_slow_clock )

!-----------------------------------------------------------------------
! average core variables on sigma surfaces
!-----------------------------------------------------------------------

  if ( dt_mean > 0 ) then
     call get_time ( Time - Time_base, secs )
     if ( mod(secs, dt_mean) == 0 ) then
        call horiz_mean (                                              &
                                                                Mgrid, &
                                                  facuv, factq, facsp, &
                                       uu(:,:,:,pres), vv(:,:,:,pres), &
                                       ta(:,:,:,pres), qv(:,:,:,pres), &
                               umean, vmean, tmean, qmean, udiv, vdiv )
     endif
  else
     umean = uu(:,:,:,pres)
     vmean = vv(:,:,:,pres)
     tmean = ta(:,:,:,pres)
     qmean = qv(:,:,:,pres)
     udiv = 0.
     vdiv = 0.
  endif
  call update_halos ( Mgrid, ps(:,:,pres) )

  call dptot ( Mgrid, dzetam, ps(:,:,pres), dmu, dmv, dmm )
  call dptot ( Mgrid, dzetaw, ps(:,:,pres), hmu, hmv, hmm )

!-----------------------------------------------------------------------
! tendency of ps
!-----------------------------------------------------------------------

  call calc_ps_tend (                                                  &
                                                                Mgrid, &
                                        ps, uu, vv, oo, dmu, dmv, dpm, &
                                                                qtend, &
                                                                 f_ps, &
                                                   past, pres, future, &
                                                                 delt )

!-----------------------------------------------------------------------
! update of ps
!-----------------------------------------------------------------------

  call step_ps (                                                       &
                                                                Mgrid, &
                                                             ps, f_ps, &
                                                         past, future, &
                                                                 delt )

!-----------------------------------------------------------------------
! diagnose dp
!-----------------------------------------------------------------------

  call dpdry ( Mgrid, dzetam, ps(:,:,future), qv(:,:,:,pres),          &
          dpu(:,:,:,future), dpv(:,:,:,future), dpm(:,:,:,future))

!-----------------------------------------------------------------------
! tendency of uu and vv
!-----------------------------------------------------------------------

  call calc_uv_tend (                                                  &
                                   Mgrid, Hmetric, Vmetric, Time_next, &
                                ps, uu, vv, oo, ta, qv, dpu, dpv, dpm, &
                                                             hpu, hpv, &
                                         dmu, dmv, dmm, hmu, hmv, hmm, &
                                             umean, vmean, udiv, vdiv, &
                                                   f_uu, f_vv, gz, fz, &
                             diffu, diffv, rdzmu, rdzwu, rdzmv, rdzwv, &
                                 cdu, cdv, dflxu, dflxv, vsatu, vsatv, &
                                             beta_u, gamma_u, ldiag_u, &
                                             beta_v, gamma_v, ldiag_v, &
               udt_mix, vdt_mix, udt_vdif, vdt_vdif, udt_adv, vdt_adv, &
                                                   udt_hpgf, vdt_hpgf, &
                                     uv_is_fluxed, past, pres, future, &
                                                                 delt )

!-----------------------------------------------------------------------
! update of uu and vv
!-----------------------------------------------------------------------

  f_uu = f_uu + udt_filt
  f_vv = f_vv + vdt_filt

  call step_uv (                                                       &
                                                     Mgrid, Time_next, &
                                                     uu, vv, dpu, dpv, &
                                                           f_uu, f_vv, &
                                                   udt_vdif, vdt_vdif, &
                                             beta_u, gamma_u, ldiag_u, &
                                             beta_v, gamma_v, ldiag_v, &
                                     uv_is_fluxed, past, pres, future, &
                                                                 delt )

!-----------------------------------------------------------------------
! diagnose sigma-dot
!-----------------------------------------------------------------------
 
  call get_divint (         Mgrid, uu(:,:,:,future), vv(:,:,:,future), &
                                          dmu, dmv, dpm(:,:,:,future), &
                                                    qtend, f_vv, f_uu )

  call get_oo ( Mgrid, zetaw, f_uu, ps(:,:,future), oo(:,:,:,future) )

  call update_halos ( Mgrid, oo(:,:,:,future) )

!-----------------------------------------------------------------------
! tendency of ta and qv
!-----------------------------------------------------------------------

  call calc_tq_tend (                                                  &
                                   Mgrid, Hmetric, Vmetric, Time_next, &
                           uu, vv, oo, ta, qv, sm, dpu, dpv, dpm, hpm, &
                                                         tmean, qmean, &
                                                     om, oh, alfomega, &
                                                           f_ta, f_qv, &
                                                               gz, fz, &
                                 diffh, rdzm, rdzw, cdh, dflxh, vsatm, &
                                                        sst, ssq, dzw, &
                                             beta_h, gamma_h, ldiag_h, &
      tdt_mix, qdt_mix, tdt_vdif, qdt_vdif, tdt_adv, qdt_adv, sdt_adv, &
                                     tq_is_fluxed, past, pres, future, &
                                                                 delt )

!-----------------------------------------------------------------------
! update of ta and qv
!-----------------------------------------------------------------------

  f_ta = f_ta + tdt_filt
  f_qv = f_qv + qdt_filt

  call step_tq (                                                       &
                                                     Mgrid, Time_next, &
                                                      ps, ta, qv, dpm, &
                                                           f_ta, f_qv, &
                                                   tdt_vdif, qdt_vdif, &
                                             beta_h, gamma_h, ldiag_h, &
                                     tq_is_fluxed, past, pres, future, &
                                                                 delt )

  call mpp_clock_end ( zetac_slow_clock )

  call mpp_clock_end ( zetac_core_driver_clock )

!-----------------------------------------------------------------------
! diagnose vertical diffusion
!-----------------------------------------------------------------------

  f_ta = gcp*dzw ! turns T_zz into theta_zz
  vbdy = sst - Tlapse*dzw(:,:,ked)
  tdt_damp = tdt_vdif  ! initially vdif + nudge
  call vdiff_eval ( Mgrid, diffh, cdh, dflxh, vsatm,                   &
                                 rdzm, rdzw, f_ta, ta, vbdy, tdt_vdif, &
                                                         pres, future )
  tdt_damp = tdt_damp - tdt_vdif
  tsurf = vbdy + sst

  f_qv = 0.
  vbdy = ssq - qlapse*dzw(:,:,ked)
  qdt_damp = qdt_vdif
  call vdiff_eval ( Mgrid, diffh, cdh, dflxh, vsatm,                   &
                                 rdzm, rdzw, f_qv, qv, vbdy, qdt_vdif, &
                                                         pres, future )
  qdt_damp = qdt_damp - qdt_vdif
  qsurf = vbdy + ssq

  f_uu = 0.
  vbdy = 0.
  udt_damp = udt_vdif
  call vdiff_eval ( Mgrid, diffu, cdu, dflxu, vsatu,                   &
                               rdzmu, rdzwu, f_uu, uu, vbdy, udt_vdif, &
                                                         pres, future )
  udt_damp = udt_damp - udt_vdif

  f_vv = 0.
  vbdy = 0.
  vdt_damp = vdt_vdif
  call vdiff_eval ( Mgrid, diffv, cdv, dflxv, vsatv,                   &
                               rdzmv, rdzwv, f_vv, vv, vbdy, vdt_vdif, &
                                                         pres, future )
  vdt_damp = vdt_damp - vdt_vdif

  call cgrid_interp_mass (Mgrid, iz,                                   &
                                  uu(:,:,:,future), vv(:,:,:,future),  &
                                                  um(:,:,:), vm(:,:,:))

!-----------------------------------------------------------------------
! diagnose kinetic energy tendency
!-----------------------------------------------------------------------

  udt_vdif = udt_vdif * dpu(:,:,:,pres)
  vdt_vdif = vdt_vdif * dpv(:,:,:,pres)

  call update_halos ( Mgrid, udt_vdif )
  call update_halos ( Mgrid, vdt_vdif )
  call cgrid_interp_mass (Mgrid, iz, udt_vdif, vdt_vdif, uatm, vatm)
  uatm = uatm / dpm(:,:,:,pres)
  vatm = vatm / dpm(:,:,:,pres)
  
  call ncdf_fms_write (Mgrid, Time_next, id_udt_vdif, gridm, uatm)
  call ncdf_fms_write (Mgrid, Time_next, id_vdt_vdif, gridm, vatm)
  call ncdf_fms_write_integral (Mgrid, Time_next, id_udt_vdif_int,     &
                                     uatm, dpm(:,:,:,pres), gfac=Grav)
  call ncdf_fms_write_integral (Mgrid, Time_next, id_vdt_vdif_int,     &
                                     vatm, dpm(:,:,:,pres), gfac=Grav)

!-----------------------------------------------------------------------
! add damping and advective diffusion to kinetic energy tendency
!-----------------------------------------------------------------------

  udt_filt = udt_filt - udt_adv
  vdt_filt = vdt_filt - vdt_adv

  udt_vdif = udt_vdif + (udt_mix + udt_damp + udt_filt) * dpu(:,:,:,pres)
  vdt_vdif = vdt_vdif + (vdt_mix + vdt_damp + vdt_filt) * dpv(:,:,:,pres)
  udt_vdif = udt_vdif * uu(:,:,:,pres)
  vdt_vdif = vdt_vdif * vv(:,:,:,pres)

  call vdifk_eval ( Mgrid, diffu, cdu, vsatu, rdzmu, rdzwu,         &
                                        uu(:,:,:,pres), uatm, utmp )
  call vdifk_eval ( Mgrid, diffv, cdv, vsatv, rdzmv, rdzwv,         &
                                        vv(:,:,:,pres), vatm, vtmp )
  udt_vdif = udt_vdif - uatm * dmu  ! correct diss heating
  vdt_vdif = vdt_vdif - vatm * dmv

  call update_halos ( Mgrid, udt_vdif )
  call update_halos ( Mgrid, vdt_vdif )
  call cgrid_interp_mass (Mgrid, iz, udt_vdif, vdt_vdif, uatm, vatm)

  kedt = (uatm + vatm) / dpm(:,:,:,pres) * (1. + qv(:,:,:,pres))
  kedt = kedt * dpm(:,:,:,pres) / dpm(:,:,:,future)

! more diagnostics

  call update_halos ( Mgrid, udt_hpgf )
  call update_halos ( Mgrid, vdt_hpgf )
  call divergence ( Mgrid, udt_hpgf, vdt_hpgf, uatm )
  call ncdf_fms_write (Mgrid, Time_next, id_divdt_hpgf, gridm, uatm)

  call cgrid_interp_mass (Mgrid, iz, udt_hpgf, vdt_hpgf, uatm, vatm)
  call ncdf_fms_write (Mgrid, Time_next, id_pgfx, gridm, uatm)
  call ncdf_fms_write (Mgrid, Time_next, id_pgfy, gridm, vatm)
  
  call ncdf_fms_write (Mgrid, Time_next, id_udt_damp, gridm, udt_damp)
  call ncdf_fms_write (Mgrid, Time_next, id_vdt_damp, gridm, vdt_damp)
  call ncdf_fms_write_integral (Mgrid, Time_next, id_udt_damp_int,     &
                                 udt_damp, dpu(:,:,:,pres), gfac=Grav)
  call ncdf_fms_write_integral (Mgrid, Time_next, id_vdt_damp_int,     &
                                 vdt_damp, dpv(:,:,:,pres), gfac=Grav)

  call get_rh ( ps(:,:,future), tsurf, qsurf, hsurf )

  call ncdf_fms_write (Mgrid, Time_next, id_tsurf, gridm, tsurf)
  call ncdf_fms_write (Mgrid, Time_next, id_qsurf, gridm, qsurf)
  call ncdf_fms_write (Mgrid, Time_next, id_hsurf, gridm, hsurf)

  return
end subroutine zetac

!#######################################################################

subroutine dpdry ( Mgrid, dz, psm, qvap, dpatu, dpatv, dpatm )

type (horiz_grid_type),                          intent (in)    ::      &
                                                                 Mgrid

real, dimension(iz),                             intent (in)    ::  dz

real, dimension(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed),              &
                                                 intent (in)    ::     &
                                                                   psm

real, dimension(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed,               &
                Mgrid%kbd:Mgrid%ked),            intent (in)    ::     &
                                                                  qvap

real, dimension(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed,               &
                Mgrid%kbd:Mgrid%ked),            intent (out)   ::     &
                                                   dpatu, dpatv, dpatm

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

  call get_dp ( Mgrid, dz, psm, dpatm )
  dpatm = dpatm/(1. + qvap)  ! dry mass
  call cgrid_interp_uv ( Mgrid, iz, dpatm, dpatu, dpatv )
  call update_halos ( Mgrid, dpatm )
  call update_halos ( Mgrid, dpatu )
  call update_halos ( Mgrid, dpatv )

end subroutine dpdry

!#######################################################################

subroutine dptot ( Mgrid, dz, psm, dmatu, dmatv, dmatm )

type (horiz_grid_type),                          intent (in)    ::     &
                                                                 Mgrid

real, dimension(iz),                             intent (in)    ::  dz

real, dimension(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed),              &
                                                 intent (in)    ::     &
                                                                   psm

real, dimension(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed,               &
                Mgrid%kbd:Mgrid%ked),            intent (out)   ::     &
                                                   dmatu, dmatv, dmatm

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

real, dimension(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed)  ::  psu, psv

! total mass

  call cgrid_interp_uv ( Mgrid, psm, psu, psv )

  call update_halos ( Mgrid, psu )
  call update_halos ( Mgrid, psv )
  call get_dp ( Mgrid, dz, psu, dmatu )
  call get_dp ( Mgrid, dz, psv, dmatv )
  call get_dp ( Mgrid, dz, psm, dmatm )

end subroutine dptot

!#######################################################################

subroutine zetac_init (                                                &
                                              Mgrid, Hmetric, Vmetric, &
                                        Time, Time_step, Time_base_in, &
                                    ppref, uuref, taref, qvref, gzref, &
                                        qc_index, nsteps, ndt_forward )

!-----------------------------------------------------------------------
! define arguments
!-----------------------------------------------------------------------

type (horiz_grid_type),                        intent (in)     ::      &
                                                                Mgrid

type (horiz_metric_type),                      intent (in)     ::      &
                                                              Hmetric

type (vert_metric_type),                       intent (in)     ::      &
                                                              Vmetric

type (time_type),                              intent (in)     ::      &
                                        Time, Time_step, Time_base_in

real, dimension(Mgrid%jbg:Mgrid%jeg,Mgrid%iz),                         &
                                                intent (in)    ::      &
                                           ppref, uuref, taref, qvref

real, dimension(Mgrid%jbg:Mgrid%jeg,Mgrid%iz),                         &
                                                intent (out)   ::      &
                                                                gzref

integer,                                        intent (in)    ::      &
                                                               nsteps

integer,                                        intent (out)   ::      &
                                                qc_index, ndt_forward

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

real, dimension(Mgrid%jbg:Mgrid%jeg) :: zfref

character (len=16) :: name
integer :: num_form, index

integer, dimension(2) :: axid_mm
integer, dimension(3) :: axid_mmm

  ibc  = Mgrid%ibc
  iec  = Mgrid%iec
  jbc  = Mgrid%jbc
  jec  = Mgrid%jec
  ibd  = Mgrid%ibd
  ied  = Mgrid%ied
  jbd  = Mgrid%jbd
  jed  = Mgrid%jed
  kbd  = Mgrid%kbd
  ked  = Mgrid%ked
  iz   = Mgrid%iz

  Time_base = Time_base_in

!-----------------------------------------------------------------------
! write version number to logfile
!-----------------------------------------------------------------------

  call write_version_number ( version, tag )

!-----------------------------------------------------------------------
! initialize clocks
!-----------------------------------------------------------------------

  zetac_core_driver_clock = mpp_clock_id( 'zetac_core'   )
  zetac_slow_clock        = mpp_clock_id( 'zetac_slow'   )
  zetac_step_clock        = mpp_clock_id( 'zetac_step'   )

!-----------------------------------------------------------------------
! allocate static atmosphere
!-----------------------------------------------------------------------

  allocate (zetam(iz), zetaw(iz))
  zetam = Vmetric%zetam
  zetaw = Vmetric%zetaw
  allocate (dzetam(iz), dzetaw(iz))
  dzetam = Vmetric%dzetam
  dzetaw = Vmetric%dzetaw

  ptop = Vmetric%ptop

!-----------------------------------------------------------------------
! allocate tendencies
!-----------------------------------------------------------------------

  allocate (f_ps(ibd:ied, jbd:jed))
  allocate (f_uu(ibd:ied, jbd:jed, kbd:ked))
  allocate (f_vv(ibd:ied, jbd:jed, kbd:ked))
  allocate (f_ta(ibd:ied, jbd:jed, kbd:ked))
  allocate (f_qv(ibd:ied, jbd:jed, kbd:ked))

  f_ps = 0.
  f_uu = 0.
  f_vv = 0.
  f_ta = 0.
  f_qv = 0.

  allocate (udt_hpgf(ibd:ied, jbd:jed, kbd:ked))
  allocate (vdt_hpgf(ibd:ied, jbd:jed, kbd:ked))

!-----------------------------------------------------------------------
! density, geopotential and pressure
!----------------------------------------- ------------------------------

  allocate (rdzm (ibd:ied, jbd:jed, kbd:ked))
  allocate (rdzw (ibd:ied, jbd:jed, kbd:ked))
  allocate (rdzmu(ibd:ied, jbd:jed, kbd:ked))
  allocate (rdzwu(ibd:ied, jbd:jed, kbd:ked))
  allocate (rdzmv(ibd:ied, jbd:jed, kbd:ked))
  allocate (rdzwv(ibd:ied, jbd:jed, kbd:ked))
  allocate (dzw  (ibd:ied, jbd:jed, kbd:ked))
  allocate (qvh  (ibd:ied, jbd:jed, kbd:ked))

!-----------------------------------------------------------------------
! allocate surface diagnostics
!-----------------------------------------------------------------------

  allocate (vbdy (ibd:ied, jbd:jed))
  allocate (dzbot(ibd:ied, jbd:jed))
  allocate (tsurf(ibd:ied, jbd:jed))
  allocate (qsurf(ibd:ied, jbd:jed))
  allocate (hsurf(ibd:ied, jbd:jed))

  qsurf = 0. ; hsurf = 0.

!-----------------------------------------------------------------------
! allocate interpolation arrays
!-----------------------------------------------------------------------

  allocate (uatm(ibd:ied, jbd:jed, kbd:ked))
  allocate (vatm(ibd:ied, jbd:jed, kbd:ked))
  allocate (utmp(ibd:ied, jbd:jed, kbd:ked))
  allocate (vtmp(ibd:ied, jbd:jed, kbd:ked))
  allocate (kpdt(ibd:ied, jbd:jed, kbd:ked))

!-----------------------------------------------------------------------
! define weights used in calculating pressure-gradient force
!-----------------------------------------------------------------------

  call press_grad_init ( Mgrid, Hmetric, Vmetric )

!-----------------------------------------------------------------------
! define coriolis and spherical or cylindrical metric parameters
!-----------------------------------------------------------------------

  call coriolis_init ( Mgrid, Hmetric )

!-----------------------------------------------------------------------
! define timestep
!-----------------------------------------------------------------------

  nsteps_leap = nsteps
  nsteps_forward = 0.5*(nsteps+1)

  call get_time ( Time_step, secs )
  delt = real(secs)

!-----------------------------------------------------------------------
! initialize forecast eq's, including open lateral boundary condition
!-----------------------------------------------------------------------

  call calc_var_tend_init (                                            &
                                        Mgrid, Hmetric, Vmetric, Time, &
                             ppref, uuref, taref, qvref, gzref, zfref, &
                                                       Tlapse, qlapse, &
                                                  facuv, factq, facsp, &
                                                                 delt, &
                                                             ndt_mean )
  
  dt_mean = max(ndt_mean,0)*secs

!-----------------------------------------------------------------------
! define tridiagonal coefficients; initialize open upper bdy condition
!-----------------------------------------------------------------------

  call step_var_init (                                                 &
                                        Mgrid, Hmetric, Vmetric, Time, &
                                           ppref, uuref, taref, qvref )


!-----------------------------------------------------------------------
! initialize vertical turbulence
!-----------------------------------------------------------------------

  call vert_turb_init ( Mgrid, Vmetric, Time, zfref, gzref, ppref )

!-----------------------------------------------------------------------
! initialize vertical diffusivity
!-----------------------------------------------------------------------

  call vert_diff_init ( Mgrid, Hmetric, Vmetric )

!-----------------------------------------------------------------------
! initialize time filter
!-----------------------------------------------------------------------

  call time_filter_init ( Mgrid, ndt_forward )

!-----------------------------------------------------------------------
! allocate for vertical diffusion
!-----------------------------------------------------------------------

  allocate (beta_u (ibc:iec, jbc:jec, kbd+1:ked))
  allocate (beta_v (ibc:iec, jbc:jec, kbd+1:ked))
  allocate (beta_h (ibc:iec, jbc:jec, kbd+1:ked))
  allocate (gamma_u(ibc:iec, jbc:jec, kbd+1:ked))
  allocate (gamma_v(ibc:iec, jbc:jec, kbd+1:ked))
  allocate (gamma_h(ibc:iec, jbc:jec, kbd+1:ked))
  allocate (ldiag_u(ibc:iec, jbc:jec, kbd+1:ked))
  allocate (ldiag_v(ibc:iec, jbc:jec, kbd+1:ked))
  allocate (ldiag_h(ibc:iec, jbc:jec, kbd+1:ked))

  gamma_u = 0.0 ; gamma_v = 0.0 ; gamma_h = 0.0

!-----------------------------------------------------------------------
! prepare flux-form transport
!-----------------------------------------------------------------------

  index = get_tracer_index ( model_atmos, prog_names(uu_) )
  call get_method ( index, 'equation', num=num_form )
  uu_is_fluxed = ( num_form == FLUX_FORM )

  index = get_tracer_index ( model_atmos, prog_names(vv_) )
  call get_method ( index, 'equation', num=num_form )
  vv_is_fluxed = ( num_form == FLUX_FORM )

  index = get_tracer_index ( model_atmos, prog_names(ta_) )
  call get_method ( index, 'equation', num=num_form )
  ta_is_fluxed = ( num_form == FLUX_FORM )

  index = get_tracer_index ( model_atmos, prog_names(qv_) )
  call get_method ( index, 'equation', num=num_form )
  qv_is_fluxed = ( num_form == FLUX_FORM )

  index = get_tracer_index ( model_atmos, prog_names(qc_) )
  call get_method ( index, 'equation', num=num_form )
  qc_is_fluxed = ( num_form == FLUX_FORM )
  qc_index = index

  uv_is_fluxed = uu_is_fluxed .or. vv_is_fluxed
  tq_is_fluxed = ta_is_fluxed .or. qv_is_fluxed

!-----------------------------------------------------------------------
! register diagnostics
!-----------------------------------------------------------------------

  axid_mmm = (/axid(xm_),axid(ym_),axid(zm_)/)
  axid_mm  = (/axid(xm_),axid(ym_)/)

  name = "pgfx"
  longname = "zonal pressure gradient force"
  units = 'm/s2'
  id_pgfx = register_diag_field (                                      &
                          model, name, axid_mmm, Time,                 &
                          longname, units, silly, range=(/-0.2, 0.2/) )

  name = "pgfy"
  longname = "meridional pressure gradient force"
  units = 'm/s2'
  id_pgfy = register_diag_field (                                      &
                          model, name, axid_mmm, Time,                 &
                          longname, units, silly, range=(/-0.2, 0.2/) )

  name = "divdt_hpgf"
  longname = "divergence tendency due to horizontal pressure gradient"
  units = 's^-2'
  id_divdt_hpgf = register_diag_field (                                &
                         model, name, axid_mmm, Time,                  &
                         longname, units, silly, range=(/-0.1, 0.1/) )

  name = "udt_vdif"
  longname = "zonal acceleration from vertical diffusion"
  units = 'm/s2'
  id_udt_vdif = register_diag_field (                                  &
                      model, name, axid_mmm, Time,                     &
                      longname, units, silly, range=(/-50., 50./) )

  name = "vdt_vdif"
  longname = "meridional acceleration from vertical diffusion"
  units = 'm/s2'
  id_vdt_vdif = register_diag_field (                                  &
                      model, name, axid_mmm, Time,                     &
                      longname, units, silly, range=(/-50., 50./) )

  name = "udt_damp"
  longname = "zonal acceleration from damping"
  units = 'm/s2'
  id_udt_damp = register_diag_field (                                  &
                      model, name, axid_mmm, Time,                     &
                      longname, units, silly, range=(/-50., 50./) )

  name = "vdt_damp"
  longname = "meridional acceleration from damping"
  units = 'm/s2'
  id_vdt_damp = register_diag_field (                                  &
                      model, name, axid_mmm, Time,                     &
                      longname, units, silly, range=(/-50., 50./) )

  name = "udt_vdif_int"
  longname = "zonal stress from vertical diffusion"
  units = 'N/m2'
  id_udt_vdif_int = register_diag_field (                              &
                      model, name, axid_mm, Time,                      &
                      longname, units, silly, range=(/-200., 200./) )

  name = "vdt_vdif_int"
  longname = "meridional stress from vertical diffusion"
  units = 'N/m2'
  id_vdt_vdif_int = register_diag_field (                              &
                      model, name, axid_mm, Time,                      &
                      longname, units, silly, range=(/-200., 200./) )

  name = "udt_damp_int"
  longname = "column zonal acceleration from damping"
  units = 'N/m2'
  id_udt_damp_int = register_diag_field (                              &
                       model, name, axid_mm, Time,                     &
                       longname, units, silly, range=(/-1.e4, 1.e4/) )

  name = "vdt_damp_int"
  longname = "column meridional acceleration from damping"
  units = 'N/m2'
  id_vdt_damp_int = register_diag_field (                              &
                       model, name, axid_mm, Time,                     &
                       longname, units, silly, range=(/-1.e4, 1.e4/) )

  name = "tsurf"
  longname = "surface temperature"
  units = 'K'
  id_tsurf = register_diag_field (                                     &
                      model, name, axid_mm, Time,                      &
                      longname, units, silly, range=(/250., 350./) )

  name = "qsurf"
  longname = "surface specific humidity"
  units = 'K'
  id_qsurf = register_diag_field (                                     &
                      model, name, axid_mm, Time,                      &
                      longname, units, silly, range=(/0., 4.0e-2/) )

  name = "hsurf"
  longname = "surface relative humidity"
  units = 'K'
  id_hsurf = register_diag_field (                                     &
                      model, name, axid_mm, Time,                      &
                      longname, units, silly, range=(/0., 1.6/) )

  return
end subroutine zetac_init

!#######################################################################

subroutine zetac_end

  return
end subroutine zetac_end

!#######################################################################

end module zetac_mod
