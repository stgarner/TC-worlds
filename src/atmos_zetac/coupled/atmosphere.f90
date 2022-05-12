module atmosphere_mod

!-----------------------------------------------------------------------
! mdodule information
!-----------------------------------------------------------------------

use fms_mod,                     only : write_version_number,          &
                                        stdout, stdlog, set_domain,    &
                                        file_exist
use mpp_mod,                     only : mpp_pe, mpp_root_pe,           &
                                        mpp_clock_id, mpp_clock_begin, &
                                        mpp_clock_end, mpp_exit,       &
                                        mpp_set_current_pelist
use mpp_io_mod,                  only : axistype, fieldtype,           &
                                        MPP_OVERWR, MPP_RDONLY,        &
                                        MPP_SINGLE, MPP_MULTI,         &
                                        MPP_NETCDF,                    &
                                        mpp_write_meta, mpp_write,     &
                                        mpp_open, mpp_close,           &
                                        mpp_get_atts, mpp_get_fields,  &
                                        mpp_get_info
use constants_mod,               only : grav, cp_air, hlv, rvgas, rdgas
use diag_manager_mod,            only : diag_axis_init,                &
                                        diag_manager_end,              &
                                        register_diag_field,           &
                                        get_base_date
use field_manager_mod,           only : model_atmos
use time_manager_mod,            only : time_type, time_manager_init,  &
                                        get_time, set_date,            &
                                        operator(+), operator(-),      &
                                        operator(.le.)
use zetac_mod,                   only : zetac, zetac_init, zetac_end
use zetac_advect_diffuse_mod,    only : advect_diffuse
use zetac_axes_mod,              only : axid, gridm
use zetac_axis_names_mod,        only : xm_, ym_, zm_
use zetac_cgrid_interp_mod,      only : cgrid_interp_mass
use zetac_check_blowup_mod,      only : check_blowup
use zetac_convert_var_mod,       only : get_rh, get_pp, get_sm
use zetac_domains_mod,           only : get_domains
use zetac_extrap_var_mod,        only : extrap_ptq, extrap_uv,         &
                                        extrap_tracer
use zetac_global_mod,            only : global
use zetac_horiz_grid_type_mod,   only : horiz_grid_type
use zetac_horiz_metric_mod,      only : horiz_metric
use zetac_horiz_metric_type_mod, only : horiz_metric_type
use zetac_history_mod,           only : write_history,                 &
                                        write_history_init
use zetac_ls_forcing_mod,        only : ls_forcing, ls_forcing_init
use zetac_moisture_driver_mod,   only : moisture_driver,               &
                                        moisture_driver_init,          &
                                        moisture_driver_end
use zetac_namelist_mod,          only : namelist_read
use zetac_ncdf_io_mod,           only : ncdf_fms_write, ncdf_io_init
use zetac_restart_mod,           only : read_restart,                  &
                                        read_restart_init,             &
                                        write_restart,                 &
                                        write_restart_init
use zetac_set_topog_mod,         only : set_topog
use zetac_step_var_mod,          only : step_tracer
use zetac_time_filter_mod,       only : time_filter
use zetac_time_pointers_mod,     only : update_time_pointers,          &
                                        get_time_pointers, ntime
use zetac_update_halos_mod,      only : update_halos, update_halos_init
use zetac_vert_metric_mod,       only : vert_metric
use zetac_vert_metric_type_mod,  only : vert_metric_type
use zetac_vert_turb_mod,         only : vert_turb

use zetac_advect_horiz_mod,      only : advect_horiz_at_mass,          &
                                        advect_horiz_quick
use zetac_advect_vert_mod,       only : advect_vert_at_mass,           &
                                        advect_vert_quick

!-----------------------------------------------------------------------
! public interfaces
!-----------------------------------------------------------------------

implicit none
private

public atmosphere_dynamics, atmosphere_state_update,                   &
       atmosphere_init, atmosphere_end

character (len=16) :: name, units
character (len=64) :: longname

integer :: atmosphere_init_clock
integer :: atmosphere_end_clock
integer :: zetac_update_clock
integer :: id_tdt_vdif, id_qdt_vdif
integer :: id_tdt_phys, id_qdt_phys
integer :: id_tdt_damp, id_qdt_damp
integer :: id_sdt_phys, id_sdt_diss, id_sdt_vdif
integer :: id_sdt_drag, id_sdt_adv, id_sdt_dry
integer :: id_tdt_vdif_int, id_qdt_vdif_int
integer :: id_tdt_phys_int, id_qdt_phys_int
integer :: id_tdt_damp_int, id_qdt_damp_int
integer :: id_sdt_phys_int, id_sdt_diss_int, id_sdt_vdif_int
integer :: id_sdt_drag_int, id_sdt_adv_int, id_sdt_dry_int
integer :: id_tdt_cond, id_qdt_cond, id_tdt_cond_int, id_qdt_cond_int
integer :: id_qdt_drag, id_qdt_drag_int
integer :: id_alfomega

!-----------------------------------------------------------------------
! version control information
!-----------------------------------------------------------------------

character(len=128)  :: version =  '$Id: atmosphere.f90,v 1.1.2.1.2.13 2005/08/07 00:12:46 stg Exp $'
character(len=128)  :: tag     =  '$Name:  $'

!-----------------------------------------------------------------------
! core namelists
!-----------------------------------------------------------------------

logical :: lxopen=.true., lyopen=.true. 
logical :: do_cartesian=.false., do_cylindrical=.false.

real :: rlonmin, rlonmax, rlatmin, rlatmax
real :: pbot, ptop

integer :: ix, iy, iz
integer :: i, j, k, ibc, iec, jbc, jec, kbd, ked
integer :: nbuf=2, nsteps=1
integer :: ndt_forward=0
integer :: qc_index

real, parameter :: silly=-99999.9
integer, parameter :: nsteps_verbose=30

real :: qvmin

!-----------------------------------------------------------------------
! static fields f(x,y)
!-----------------------------------------------------------------------

real, allocatable, dimension(:,:) :: topog, sst, ssq

!-----------------------------------------------------------------------
! static fields f(y,z) (2D default for external fields)
!-----------------------------------------------------------------------

real, allocatable, dimension(:,:) :: ppref
real, allocatable, dimension(:,:) :: uuref
real, allocatable, dimension(:,:) :: taref
real, allocatable, dimension(:,:) :: qvref
real, allocatable, dimension(:,:) :: gzref

!-----------------------------------------------------------------------
! prognostic variables
!-----------------------------------------------------------------------

real, allocatable, dimension(:,:,:)   :: ps
real, allocatable, dimension(:,:,:,:) :: uu
real, allocatable, dimension(:,:,:,:) :: vv
real, allocatable, dimension(:,:,:,:) :: oo
real, allocatable, dimension(:,:,:,:) :: ta
real, allocatable, dimension(:,:,:,:) :: qv
real, allocatable, dimension(:,:,:,:) :: qc
real, allocatable, dimension(:,:,:,:) :: sm
real, allocatable, dimension(:,:,:,:) :: dpu, dpv, dpm
real, allocatable, dimension(:,:,:) :: hpu, hpv, hpm
real, allocatable, dimension(:,:,:) :: dmu, dmv, dmm
real, allocatable, dimension(:,:,:) :: hmu, hmv, hmm

!-----------------------------------------------------------------------
! interpolated variables
!----------------------------------------- ------------------------------

real, allocatable, dimension(:,:,:) :: om, oh, um, vm, pf, ph
real, allocatable, dimension(:,:,:) :: alfomega

!-----------------------------------------------------------------------
! tracer variables
!-----------------------------------------------------------------------

real, allocatable, dimension(:,:,:) :: f_qc
real, allocatable, dimension(:,:,:) :: uatm, vatm

!-----------------------------------------------------------------------
! diagnostic variables
!-----------------------------------------------------------------------

real, allocatable, dimension(:,:)   :: rain
real, allocatable, dimension(:,:)   :: snow
real, allocatable, dimension(:,:,:) :: gz, fz, rh, th, qh

!-----------------------------------------------------------------------
! diffusion variables
!-----------------------------------------------------------------------

real, allocatable, dimension(:,:,:) :: diffu, diffv, diffh
real, allocatable, dimension(:,:)   :: dflxu, dflxv, dflxh
real, allocatable, dimension(:,:)   :: cdu, cdv, cdh
real, allocatable, dimension(:,:)   :: vsatu, vsatv, vsatm

!-----------------------------------------------------------------------
! tendency variables
!-----------------------------------------------------------------------

real, allocatable, dimension(:,:)   :: tdt, qdt, sdt
real, allocatable, dimension(:,:,:) :: ttend, qtend
real, allocatable, dimension(:,:,:) :: kedt, pedt, zqdt
real, allocatable, dimension(:,:,:) :: tdt_phys, qdt_phys
real, allocatable, dimension(:,:,:) :: sdt_phys, sdt_diss, sdt_vdif
real, allocatable, dimension(:,:,:) :: sdt_drag, sdt_dry
real, allocatable, dimension(:,:,:) :: tdt_cond, qdt_cond
real, allocatable, dimension(:,:,:) :: udt_damp, vdt_damp
real, allocatable, dimension(:,:,:) :: tdt_damp, qdt_damp
real, allocatable, dimension(:,:,:) :: udt_vdif, vdt_vdif
real, allocatable, dimension(:,:,:) :: tdt_vdif, qdt_vdif
real, allocatable, dimension(:,:,:) :: udt_filt, vdt_filt
real, allocatable, dimension(:,:,:) :: tdt_filt, qdt_filt
real, allocatable, dimension(:,:,:) :: udt_mix, vdt_mix
real, allocatable, dimension(:,:,:) :: tdt_mix, qdt_mix
real, allocatable, dimension(:,:,:) :: udt_adv, vdt_adv
real, allocatable, dimension(:,:,:) :: tdt_adv, qdt_adv, sdt_adv, pdt_adv
real, allocatable, dimension(:,:,:) :: umean, vmean, tmean, qmean
real, allocatable, dimension(:,:,:) :: udiv, vdiv, var

type (horiz_grid_type),   save :: Mgrid
type (horiz_metric_type), save :: Hmetric
type (vert_metric_type),  save :: Vmetric

type (time_type),         save :: Time_base
type (time_type),         save :: Time_prev
type (time_type),         save :: Time_next
type (time_type),         save :: Time_init
type (time_type),         save :: Time_step
type (time_type),         save :: Time_last

!-----------------------------------------------------------------------
! metadata variables
!-----------------------------------------------------------------------

character (len=16), parameter :: model='atmos_mod'

!-----------------------------------------------------------------------
! timeipointers and controls
!-----------------------------------------------------------------------

logical :: early_exit, do_forward=.false., do_coldstart=.false.
integer :: past, pres, future
integer :: secs, secs_forward, na=0
real    :: run_time, elapsed_time, percent, mxval=0.0, delt
real    :: srad, sdif, sdis, sadv, sdrg, sdry, stot

contains

!#######################################################################

subroutine atmosphere_dynamics ( Time, Time_step )

type (time_type),                               intent(in)     ::      &
                                                      Time, Time_step

integer :: i, j, k

!-----------------------------------------------------------------------
! reset time-types -- use half timestep for forward step
!-----------------------------------------------------------------------

  if ( do_coldstart .or. do_forward ) then
     Time_prev = Time
     do_coldstart = .false.
  else
     Time_prev = Time - Time_step
  endif

  Time_next = Time + Time_step

!-----------------------------------------------------------------------
! get current time, timestep and time indices
!-----------------------------------------------------------------------

  call get_time ( Time_next - Time_prev, secs )
  delt = real(secs)

  call get_time_pointers ( past, pres, future )

!-----------------------------------------------------------------------
! vertical diffusivity
!-----------------------------------------------------------------------

  call vert_turb (                                                     &
                                                          Mgrid, Time, &
                                   uu(:,:,ked,pres), vv(:,:,ked,pres), &
                     ps(:,:,pres), ta(:,:,ked,pres), qv(:,:,ked,pres), &
                                                  diffu, diffv, diffh, &
              cdu, cdv, cdh, dflxu, dflxv, dflxh, vsatu, vsatv, vsatm, &
                                                                 delt )

!-----------------------------------------------------------------------
! KE tendency from (lagged) qv tendency
!-----------------------------------------------------------------------

  udt_mix = 0.5*uu(:,:,:,pres)**2*dpu(:,:,:,pres)
  vdt_mix = 0.5*vv(:,:,:,pres)**2*dpv(:,:,:,pres)
  call cgrid_interp_mass ( Mgrid, iz, udt_mix, vdt_mix, uatm, vatm )
  zqdt = (uatm + vatm) / dpm(:,:,:,pres) * qtend

!-----------------------------------------------------------------------
! dynamics update of core variables
!-----------------------------------------------------------------------

  call zetac (                                                         &
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
! advection-diffusion update of tracer
!-----------------------------------------------------------------------

  call advect_diffuse (                                                &
                                              Mgrid, Hmetric, Vmetric, &
                                   uu, vv, oo, qc, dpu, dpv, dpm, hpm, &
                                                                 f_qc, &
                                         qc_index, past, pres, future, &
                                                                 delt )

  call step_tracer (                                                   &
                                                                Mgrid, &
                                                             qc, f_qc, &
                                                         past, future, &
                                                                 delt )
  qc(:,:,:,future) = max(0., qc(:,:,:,future))

! update ttend, qtend with time filter alteration and external forcing

  ttend = tdt_vdif + tdt_damp + tdt_mix + tdt_filt
  qtend = qdt_vdif + qdt_damp + qdt_mix + qdt_filt

  tdt_damp = tdt_damp + tdt_filt  ! add time filter to damping
  qdt_damp = qdt_damp + qdt_filt

  return
end subroutine atmosphere_dynamics

!#######################################################################

subroutine atmosphere_state_update ( Time )

type (time_type),  intent(in) :: Time
integer :: i, j

  call mpp_clock_begin ( zetac_update_clock )

!-----------------------------------------------------------------------
! large-scale forcing
!-----------------------------------------------------------------------

  call update_halos ( Mgrid, ttend )
  call update_halos ( Mgrid, qtend )

  tdt_phys = ta(:,:,:,future)
  qdt_phys = qv(:,:,:,future)

  call ls_forcing (                                                    &
                                               Mgrid, Time, Time_next, &
                                   ta(:,:,:,future), qv(:,:,:,future), &
                                                                 delt )

  tdt_phys = (ta(:,:,:,future) - tdt_phys)/delt
  qdt_phys = (qv(:,:,:,future) - qdt_phys)/delt

!-----------------------------------------------------------------------
! microphysics and cloud scheme
!-----------------------------------------------------------------------

  tdt_cond = ta(:,:,:,future)
  qdt_cond = qv(:,:,:,future)

  call moisture_driver (                                               &
                                              Mgrid, Hmetric, Vmetric, &
                                                                 Time, &
                                     ps(:,:,future), ta(:,:,:,future), &
                                   qv(:,:,:,future), qc(:,:,:,future), &
                                                           rain, snow, &
                                                                 delt )

  call update_halos ( Mgrid, ta(:,:,:,future) )

!-----------------------------------------------------------------------
! extrapolate outside open boundaries
!-----------------------------------------------------------------------

  call extrap_ptq (                                                    &
                                                                Mgrid, &
                   ps(:,:,future), ta(:,:,:,future), qv(:,:,:,future), &
                                           ppref(:,ked), taref, qvref, &
				                                 delt )

  call extrap_uv (                                                     &
                                                                Mgrid, &
                                   uu(:,:,:,future), vv(:,:,:,future), &
                                                                uuref )

  call extrap_tracer (                                                 &
                                                                Mgrid, &
                                                     qc(:,:,:,future), &
                                                                 delt )

!-----------------------------------------------------------------------
! apply halo updates to core variables
!-----------------------------------------------------------------------

  call update_halos (                                                  &
                                                                Mgrid, &
                                               ps, uu, vv, ta, qv, qc, &
                                                               future )

  tdt_cond = (ta(:,:,:,future) - tdt_cond)/delt
  qdt_cond = (qv(:,:,:,future) - qdt_cond)/delt

!-----------------------------------------------------------------------
! update ttend, qtend with radiation and microphysics
!-----------------------------------------------------------------------

  ttend = ttend + tdt_phys + tdt_cond
  qtend = qtend + qdt_phys + qdt_cond

!  tdt_phys = Cp_air*(tdt_phys + tdt_damp + tdt_filt)
!  qdt_phys =    Hlv*(qdt_phys + qdt_damp + qdt_filt)
  tdt_phys = Cp_air*(tdt_phys + tdt_damp)
  qdt_phys =    Hlv*(qdt_phys + qdt_damp)
  tdt_vdif = Cp_air*(tdt_vdif + tdt_mix)
  qdt_vdif =    Hlv*(qdt_vdif + qdt_mix)

!-----------------------------------------------------------------------
! energy tendencies
!-----------------------------------------------------------------------

  pedt = (ttend*Cp_air + qtend*Hlv) * dpm(:,:,:,pres)/                 &
                                      dpm(:,:,:,future)
  
  zqdt = (zqdt + fz*qtend) * dpm(:,:,:,pres)/                          &
                             dpm(:,:,:,future)    ! precipitation drag
  call update_halos ( Mgrid, zqdt )

  kedt = kedt + zqdt
  pedt = pedt - kedt
  ta(:,:,:,future) = ta(:,:,:,future) - kedt/Cp_air*delt

  call update_halos ( Mgrid, ta(:,:,:,future) )

!-----------------------------------------------------------------------
! time-filter tendencies and qv limiter
!-----------------------------------------------------------------------

  call time_filter (                                                   &
                                                                Mgrid, &
                                                       uu, vv, ta, qv, &
                               udt_filt, vdt_filt, tdt_filt, qdt_filt, &
                                                   past, pres, future )

  do k=2,ked-2
     cdu = dim(qvmin,qv(:,:,k,pres) + qdt_filt(:,:,k))
     qdt_filt(:,:,k  ) = qdt_filt(:,:,k  ) + cdu
     qdt_filt(:,:,k+1) = qdt_filt(:,:,k+1) - cdu*                      &
            dpm(:,:,k,pres)/dpm(:,:,k+1,pres)
  enddo

  udt_filt = udt_filt/delt
  vdt_filt = vdt_filt/delt
  tdt_filt = tdt_filt/delt
  qdt_filt = qdt_filt/delt

  call get_pp ( ps(:,:,future), Vmetric%zetam, ptop, pf )
  call get_rh ( pf, ta(:,:,:,future), qv(:,:,:,future), rh, 1.e-6 )
  call get_sm (ta(:,:,:,future), pf, qv(:,:,:,future), rh, sm(:,:,:,future))

! entropy tendencies
  
  sdt_phys = (tdt_phys + qdt_phys)/ta(:,:,:,pres)  ! radiation & damping
  sdt_vdif = (tdt_vdif + qdt_vdif)/ta(:,:,:,pres)  ! vertical diffusion
  sdt_diss = -kedt/ta(:,:,:,pres)                  ! dissipative heating
  sdt_drag = -zqdt/ta(:,:,:,pres)                  ! precip drag

! run-time diagnostics

  if (mod(na,nsteps_verbose) == 0) then

  call global (Mgrid, sdt_drag(ibc:iec,jbc:jec,kbd+1:ked), sdrg, &
                           dpm(ibc:iec,jbc:jec,kbd+1:ked,pres))
  if (mpp_pe() == 0) write (6,'("sdrg", 3f10.3)') sdrg, sdry, sadv+sdry

  call global (Mgrid, sdt_phys(ibc:iec,jbc:jec,kbd+1:ked), srad, &
                           dpm(ibc:iec,jbc:jec,kbd+1:ked,pres))
  call global (Mgrid, sdt_vdif(ibc:iec,jbc:jec,kbd+1:ked), sdif, &
                           dpm(ibc:iec,jbc:jec,kbd+1:ked,pres))
  call global (Mgrid, sdt_diss(ibc:iec,jbc:jec,kbd+1:ked), sdis, &
                           dpm(ibc:iec,jbc:jec,kbd+1:ked,pres))

  if (mpp_pe() == 0) write (6,'("stot", 4f10.3)') srad, sdif, sdis, stot

  if (mpp_pe() == 0) write (6,'("entropy error", 2f10.3)') &
                              sadv+sdry-stot, (sadv+sdry-stot)/abs(srad)

endif

!-----------------------------------------------------------------------
! is the next timestep a forward step?
!-----------------------------------------------------------------------

  if ( ndt_forward > 0 ) then
     call get_time ( Time_next - Time_base, secs )
     do_forward = ( mod(secs,secs_forward) == 0 )
  endif

!-----------------------------------------------------------------------
! permute time indices
!-----------------------------------------------------------------------

  call update_time_pointers ( past, pres, future, do_forward )

  call mpp_clock_end ( zetac_update_clock )

!-----------------------------------------------------------------------
! is solution blowing up?
!-----------------------------------------------------------------------

  call check_blowup ( Mgrid, oo(:,:,:,pres), ta(:,:,:,pres),           &
                                                    mxval, early_exit )

  call get_time ( Time_next - Time_init, secs )
  elapsed_time = real(secs)

  percent = 100.0*elapsed_time/run_time
  early_exit = ( early_exit .and. percent < 100.0 )

!-----------------------------------------------------------------------
! soft crash
!-----------------------------------------------------------------------

  if ( early_exit ) then

     write (stdout(),                                                  &
           '(/,"exiting early with test variable beyond tolerance",/)')
     mxval = 0.

     do while ( Time_next .le. Time_last )
        call get_time ( Time_next - Time_init, secs )
        elapsed_time = real(secs)
        percent = 100.0*elapsed_time/run_time
        na=na+1
!        write (stdout(), '("na =",i6,6x,"time =",f10.0,2x,             &
!&       "(",f6.1,"% )",6x,"maxval = *****")') na, elapsed_time, percent
        call write_history (                                           &
                                 Mgrid, Hmetric, Time_next, Time_step, &
      ps(:,:,pres), dpu(:,:,:,pres), dpv(:,:,:,pres), dpm(:,:,:,pres), &
                                       uu(:,:,:,pres), vv(:,:,:,pres), &
                                       oo(:,:,:,pres), ta(:,:,:,pres), &
                       qv(:,:,:,pres), qc(:,:,:,pres), sm(:,:,:,pres), &
                                             om, um, vm, ttend, qtend, &
                                               rh, gz, fz, kedt, pedt, &
                                    ppref, uuref, taref, qvref, gzref, &
				         elapsed_time, percent, mxval, &
                                                                 delt )
        Time_next = Time_next + Time_step
        call ncdf_fms_write (Mgrid, Time_next, id_tdt_damp, gridm, tdt_damp)
        call ncdf_fms_write (Mgrid, Time_next, id_qdt_damp, gridm, qdt_damp)
        call ncdf_fms_write (Mgrid, Time_next, id_tdt_vdif, gridm, tdt_vdif)
        call ncdf_fms_write (Mgrid, Time_next, id_qdt_vdif, gridm, qdt_vdif)
        call ncdf_fms_write (Mgrid, Time_next, id_tdt_phys, gridm, tdt_phys)
        call ncdf_fms_write (Mgrid, Time_next, id_qdt_phys, gridm, qdt_phys)
        call ncdf_fms_write (Mgrid, Time_next, id_sdt_phys, gridm, sdt_phys)
        call ncdf_fms_write (Mgrid, Time_next, id_sdt_vdif, gridm, sdt_vdif)
        call ncdf_fms_write (Mgrid, Time_next, id_sdt_diss, gridm, sdt_diss)
        call ncdf_fms_write (Mgrid, Time_next, id_sdt_adv, gridm, sdt_adv)
        call ncdf_fms_write (Mgrid, Time_next, id_sdt_dry, gridm, sdt_dry)
        call ncdf_fms_write (Mgrid, Time_next, id_sdt_drag, gridm, sdt_drag)
        call ncdf_fms_write (Mgrid, Time_next, id_tdt_cond, gridm, tdt_cond)
        call ncdf_fms_write (Mgrid, Time_next, id_qdt_cond, gridm, qdt_cond)
        call ncdf_fms_write (Mgrid, Time_next, id_qdt_drag, gridm, zqdt)
        call ncdf_fms_write (Mgrid, Time_next, id_alfomega, gridm, alfomega)
     enddo

     call diag_manager_end  ( Time_last )
     call atmosphere_end    ( Time_last )
     call mpp_exit
     stop

  else

!-----------------------------------------------------------------------
! send physics data to diagnostics
!-----------------------------------------------------------------------

     na=na+1
!     write (stdout(), '("na =",i6,6x,"time =",f10.0,2x,"(",f6.1,"% )", &
!&               6x,"maxval =",f6.2)')  na, elapsed_time, percent, mxval

!-----------------------------------------------------------------------
! send dycore data to diagnostics
!-----------------------------------------------------------------------

     call write_history (                                              &
                                 Mgrid, Hmetric, Time_next, Time_step, &
      ps(:,:,pres), dpu(:,:,:,pres), dpv(:,:,:,pres), dpm(:,:,:,pres), &
                                       uu(:,:,:,pres), vv(:,:,:,pres), &
                                       oo(:,:,:,pres), ta(:,:,:,pres), &
                       qv(:,:,:,pres), qc(:,:,:,pres), sm(:,:,:,pres), &
                                             om, um, vm, ttend, qtend, &
                                               rh, gz, fz, kedt, pedt, &
                                    ppref, uuref, taref, qvref, gzref, &
                                         elapsed_time, percent, mxval, &
                                                                 delt )

!-----------------------------------------------------------------------
! send physics data to diagnostics
!-----------------------------------------------------------------------

     call ncdf_fms_write (Mgrid, Time_next, id_tdt_damp, gridm, tdt_damp)
     call ncdf_fms_write (Mgrid, Time_next, id_qdt_damp, gridm, qdt_damp)
     call ncdf_fms_write (Mgrid, Time_next, id_tdt_vdif, gridm, tdt_vdif)
     call ncdf_fms_write (Mgrid, Time_next, id_qdt_vdif, gridm, qdt_vdif)
     call ncdf_fms_write (Mgrid, Time_next, id_tdt_phys, gridm, tdt_phys)
     call ncdf_fms_write (Mgrid, Time_next, id_qdt_phys, gridm, qdt_phys)
     call ncdf_fms_write (Mgrid, Time_next, id_sdt_phys, gridm, sdt_phys)
     call ncdf_fms_write (Mgrid, Time_next, id_sdt_vdif, gridm, sdt_vdif)
     call ncdf_fms_write (Mgrid, Time_next, id_sdt_diss, gridm, sdt_diss)
     call ncdf_fms_write (Mgrid, Time_next, id_sdt_adv, gridm, sdt_adv)
     call ncdf_fms_write (Mgrid, Time_next, id_sdt_dry, gridm, sdt_dry)
     call ncdf_fms_write (Mgrid, Time_next, id_sdt_drag, gridm, sdt_drag)
     call ncdf_fms_write (Mgrid, Time_next, id_tdt_cond, gridm, tdt_cond)
     call ncdf_fms_write (Mgrid, Time_next, id_qdt_cond, gridm, qdt_cond)
     call ncdf_fms_write (Mgrid, Time_next, id_qdt_drag, gridm, zqdt)
     call ncdf_fms_write (Mgrid, Time_next, id_alfomega, gridm, alfomega)

     do j=jbc,jec
        do i=ibc,iec
           tdt(i,j) = Cp_air/Grav * sum (tdt_damp(i,j,kbd+1:ked)*      &
                                              dpm(i,j,kbd+1:ked,pres))
           qdt(i,j) = Hlv/Grav * sum (qdt_damp(i,j,kbd+1:ked)*         &
                                           dpm(i,j,kbd+1:ked,pres))
           sdt(i,j) = sum (sdt_drag(i,j,kbd+1:ked)*                    &
                                dpm(i,j,kbd+1:ked,pres)) / Grav
        enddo
     enddo

     call ncdf_fms_write (Mgrid, Time_next, id_tdt_damp_int, gridm, tdt)
     call ncdf_fms_write (Mgrid, Time_next, id_qdt_damp_int, gridm, qdt)
     call ncdf_fms_write (Mgrid, Time_next, id_sdt_drag_int, gridm, sdt)

     do j=jbc,jec
        do i=ibc,iec
           tdt(i,j) = sum (tdt_phys(i,j,kbd+1:ked)*                    &
                                dpm(i,j,kbd+1:ked,pres)) / Grav
           qdt(i,j) = sum (qdt_phys(i,j,kbd+1:ked)*                    &
                                dpm(i,j,kbd+1:ked,pres)) / Grav
           sdt(i,j) = sum (sdt_phys(i,j,kbd+1:ked)*                    &
                                dpm(i,j,kbd+1:ked,pres)) / Grav
        enddo
     enddo

     call ncdf_fms_write (Mgrid, Time_next, id_tdt_phys_int, gridm, tdt)
     call ncdf_fms_write (Mgrid, Time_next, id_qdt_phys_int, gridm, qdt)
     call ncdf_fms_write (Mgrid, Time_next, id_sdt_phys_int, gridm, sdt)

     do j=jbc,jec
        do i=ibc,iec
           tdt(i,j) = Cp_air/Grav * sum (tdt_cond(i,j,kbd+1:ked)*      &
                                              dpm(i,j,kbd+1:ked,pres))
           qdt(i,j) = Hlv/Grav * sum (qdt_cond(i,j,kbd+1:ked)*         &
                                           dpm(i,j,kbd+1:ked,pres))
           sdt(i,j) = sum (sdt_diss(i,j,kbd+1:ked)*                    &
                                dpm(i,j,kbd+1:ked,pres)) / Grav
        enddo
     enddo

     call ncdf_fms_write (Mgrid, Time_next, id_tdt_cond_int, gridm, tdt)
     call ncdf_fms_write (Mgrid, Time_next, id_qdt_cond_int, gridm, qdt)
     call ncdf_fms_write (Mgrid, Time_next, id_sdt_diss_int, gridm, sdt)

     do j=jbc,jec
        do i=ibc,iec
           tdt(i,j) = sum (tdt_vdif(i,j,kbd+1:ked)*                    &
                                dpm(i,j,kbd+1:ked,pres)) / Grav
           qdt(i,j) = sum (qdt_vdif(i,j,kbd+1:ked)*                    &
                                dpm(i,j,kbd+1:ked,pres)) / Grav
           sdt(i,j) = sum (sdt_vdif(i,j,kbd+1:ked)*                    &
                                dpm(i,j,kbd+1:ked,pres)) / Grav
        enddo
     enddo
     
     call ncdf_fms_write (Mgrid, Time_next, id_tdt_vdif_int, gridm, tdt)
     call ncdf_fms_write (Mgrid, Time_next, id_qdt_vdif_int, gridm, qdt)
     call ncdf_fms_write (Mgrid, Time_next, id_sdt_vdif_int, gridm, sdt)

     do j=jbc,jec
        do i=ibc,iec
           qdt(i,j) = sum (zqdt(i,j,kbd+1:ked)*                        &
                            dpm(i,j,kbd+1:ked,pres)) / Grav
           sdt(i,j) = sum (sdt_adv(i,j,kbd+1:ked)                      &
                              *dpm(i,j,kbd+1:ked,pres)) / Grav
        enddo
     enddo

     call ncdf_fms_write (Mgrid, Time_next, id_qdt_drag_int, gridm, qdt)
     call ncdf_fms_write (Mgrid, Time_next, id_sdt_adv_int, gridm, sdt)

     do j=jbc,jec
        do i=ibc,iec
           sdt(i,j) = sum (sdt_dry(i,j,kbd+1:ked)*                     &
                               dpm(i,j,kbd+1:ked,pres)) / Grav
        enddo
     enddo

     call ncdf_fms_write (Mgrid, Time_next, id_sdt_dry_int, gridm, sdt)

  endif

  return
end subroutine atmosphere_state_update

!#######################################################################

subroutine atmosphere_init (                                           &
                    Time, Time_init_arg, Time_step_arg, Time_last_arg )

type (time_type),            intent (in)                       ::      & 
                    Time, Time_init_arg, Time_step_arg, Time_last_arg

integer :: date(6)
integer :: axid_mm(2), axid_mmm(3)
real, dimension(20) :: zeta
real, dimension(:), allocatable :: zetaw

!-----------------------------------------------------------------------
! write version number to logfile
!-----------------------------------------------------------------------

  call write_version_number (version, tag)

!-----------------------------------------------------------------------
! initialize clocks
!-----------------------------------------------------------------------

  atmosphere_init_clock       = mpp_clock_id( 'atmosphere_init'       )
  atmosphere_end_clock        = mpp_clock_id( 'atmosphere_end'        )
  zetac_update_clock          = mpp_clock_id( 'zetac_atmos_update'    )

  call mpp_clock_begin ( atmosphere_init_clock )

!-----------------------------------------------------------------------
! store time_types
!-----------------------------------------------------------------------

  Time_init = Time_init_arg
  Time_step = Time_step_arg
  Time_last = Time_last_arg

!-----------------------------------------------------------------------
! read 'grid' and 'region' namelists
!-----------------------------------------------------------------------

  call namelist_read (                                                 &
                                                           ix, iy, iz, &
                                                         nsteps, nbuf, &
                         lxopen, lyopen, do_cartesian, do_cylindrical, &
                                   rlonmin, rlonmax, rlatmin, rlatmax, &
                                                     pbot, ptop, zeta )

   allocate (zetaw(iz))
   zetaw = zeta(1:iz)

!-----------------------------------------------------------------------
! do domain decomposition and define grid_type variable
!-----------------------------------------------------------------------

  call get_domains (                                                   &
                                                     ix, iy, iz, nbuf, &
                                                Mgrid, lxopen, lyopen )

  ibc = Mgrid%ibc
  iec = Mgrid%iec
  jbc = Mgrid%jbc
  jec = Mgrid%jec
  kbd = Mgrid%kbd
  ked = Mgrid%ked

!-----------------------------------------------------------------------
! initialize wrapper for mpp_update_domains
!-----------------------------------------------------------------------

  call update_halos_init ( Mgrid )

!-----------------------------------------------------------------------
! allocate static variable arrays
!-----------------------------------------------------------------------

  allocate (topog(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed))
  allocate (  sst(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed))
  allocate (  ssq(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed))

!-----------------------------------------------------------------------
! allocate reference fields
!-----------------------------------------------------------------------

  allocate (ppref(Mgrid%jbg:Mgrid%jeg, Mgrid%kbd:Mgrid%ked))
  allocate (uuref(Mgrid%jbg:Mgrid%jeg, Mgrid%kbd:Mgrid%ked))
  allocate (taref(Mgrid%jbg:Mgrid%jeg, Mgrid%kbd:Mgrid%ked))
  allocate (qvref(Mgrid%jbg:Mgrid%jeg, Mgrid%kbd:Mgrid%ked))
  allocate (gzref(Mgrid%jbg:Mgrid%jeg, Mgrid%kbd:Mgrid%ked))

!-----------------------------------------------------------------------
! allocate prognostic arrays
!-----------------------------------------------------------------------

  allocate (ps(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,               &
                                    ntime))
  allocate (uu(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,               &
                                    Mgrid%kbd:Mgrid%ked, ntime))
  allocate (vv(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,               &
                                    Mgrid%kbd:Mgrid%ked, ntime))
  allocate (oo(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,               &
                                    Mgrid%kbd:Mgrid%ked, ntime))
  allocate (ta(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,               &
                                    Mgrid%kbd:Mgrid%ked, ntime))
  allocate (qv(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,               &
                                    Mgrid%kbd:Mgrid%ked, ntime))
  allocate (qc(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,               &
                                    Mgrid%kbd:Mgrid%ked, ntime))
  allocate (sm(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,               &
                                    Mgrid%kbd:Mgrid%ked, ntime))

!-----------------------------------------------------------------------
! allocate density fields
!-----------------------------------------------------------------------

  allocate (dpu(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                                     Mgrid%kbd:Mgrid%ked, ntime))
  allocate (dpv(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                                     Mgrid%kbd:Mgrid%ked, ntime))
  allocate (dpm(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                                     Mgrid%kbd:Mgrid%ked, ntime))
  allocate (hpu(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                                     Mgrid%kbd:Mgrid%ked))
  allocate (hpv(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                                     Mgrid%kbd:Mgrid%ked))
  allocate (hpm(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                                     Mgrid%kbd:Mgrid%ked))
  allocate (dmu(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                                     Mgrid%kbd:Mgrid%ked))
  allocate (dmv(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                                     Mgrid%kbd:Mgrid%ked))
  allocate (dmm(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                                     Mgrid%kbd:Mgrid%ked))
  allocate (hmu(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                                     Mgrid%kbd:Mgrid%ked))
  allocate (hmv(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                                     Mgrid%kbd:Mgrid%ked))
  allocate (hmm(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                                     Mgrid%kbd:Mgrid%ked))

!-----------------------------------------------------------------------
! allocate interpolated variables
!-----------------------------------------------------------------------

  allocate (om(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,               &
                                    Mgrid%kbd:Mgrid%ked))
  allocate (oh(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,               &
                                    Mgrid%kbd:Mgrid%ked))
  allocate (um(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,               &
                                    Mgrid%kbd:Mgrid%ked))
  allocate (vm(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,               &
                                    Mgrid%kbd:Mgrid%ked))
  allocate (pf(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,               &
                                    Mgrid%kbd:Mgrid%ked))
  allocate (ph(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,               &
                                    Mgrid%kbd:Mgrid%ked))
  allocate (gz(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,               &
                                    Mgrid%kbd:Mgrid%ked))
  allocate (fz(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,               &
                                    Mgrid%kbd:Mgrid%ked))
  allocate (rh(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,               &
                                    Mgrid%kbd:Mgrid%ked))
  allocate (th(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,               &
                                    Mgrid%kbd:Mgrid%ked))
  allocate (qh(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,               &
                                    Mgrid%kbd:Mgrid%ked))
  allocate (alfomega(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,         &
                                          Mgrid%kbd:Mgrid%ked))

!-----------------------------------------------------------------------
! allocate tendency arrays
!-----------------------------------------------------------------------

  allocate (ttend(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,            &
                                       Mgrid%kbd:Mgrid%ked))
  allocate (qtend(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,            &
                                       Mgrid%kbd:Mgrid%ked))
  allocate (tdt  (Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed))
  allocate (qdt  (Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed))
  allocate (sdt  (Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed))
  allocate (tdt_phys(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,         &
                                          Mgrid%kbd:Mgrid%ked))
  allocate (qdt_phys(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,         &
                                          Mgrid%kbd:Mgrid%ked))
  allocate (sdt_phys(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,         &
                                          Mgrid%kbd:Mgrid%ked))
  allocate (sdt_vdif(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,         &
                                          Mgrid%kbd:Mgrid%ked))
  allocate (sdt_diss(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,         &
                                          Mgrid%kbd:Mgrid%ked))
  allocate (sdt_drag(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,         &
                                          Mgrid%kbd:Mgrid%ked))
  allocate (tdt_cond(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,         &
                                          Mgrid%kbd:Mgrid%ked))
  allocate (qdt_cond(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,         &
                                          Mgrid%kbd:Mgrid%ked))
  allocate (udt_damp(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,         &
                                          Mgrid%kbd:Mgrid%ked))
  allocate (vdt_damp(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,         &
                                          Mgrid%kbd:Mgrid%ked))
  allocate (tdt_damp(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,         &
                                          Mgrid%kbd:Mgrid%ked))
  allocate (qdt_damp(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,         &
                                          Mgrid%kbd:Mgrid%ked))
  allocate (udt_vdif(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,         &
                                          Mgrid%kbd:Mgrid%ked))
  allocate (vdt_vdif(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,         &
                                          Mgrid%kbd:Mgrid%ked))
  allocate (tdt_vdif(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,         &
                                          Mgrid%kbd:Mgrid%ked))
  allocate (qdt_vdif(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,         &
                                          Mgrid%kbd:Mgrid%ked))
  allocate (udt_filt(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,         &
                                          Mgrid%kbd:Mgrid%ked))
  allocate (vdt_filt(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,         &
                                          Mgrid%kbd:Mgrid%ked))
  allocate (tdt_filt(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,         &
                                          Mgrid%kbd:Mgrid%ked))
  allocate (qdt_filt(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,         &
                                          Mgrid%kbd:Mgrid%ked))
  allocate (udt_mix(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,          &
                                         Mgrid%kbd:Mgrid%ked))
  allocate (vdt_mix(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,          &
                                         Mgrid%kbd:Mgrid%ked))
  allocate (tdt_mix(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,          &
                                         Mgrid%kbd:Mgrid%ked))
  allocate (qdt_mix(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,          &
                                         Mgrid%kbd:Mgrid%ked))
  allocate (udt_adv(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,          &
                                         Mgrid%kbd:Mgrid%ked))
  allocate (vdt_adv(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,          &
                                         Mgrid%kbd:Mgrid%ked))
  allocate (tdt_adv(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,          &
                                         Mgrid%kbd:Mgrid%ked))
  allocate (qdt_adv(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,          &
                                         Mgrid%kbd:Mgrid%ked))
  allocate (sdt_adv(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,          &
                                         Mgrid%kbd:Mgrid%ked))
  allocate (pdt_adv(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,          &
                                         Mgrid%kbd:Mgrid%ked))
  allocate (sdt_dry(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,          &
                                         Mgrid%kbd:Mgrid%ked))
  allocate (kedt(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,             &
                                      Mgrid%kbd:Mgrid%ked))
  allocate (pedt(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,             &
                                      Mgrid%kbd:Mgrid%ked))
  allocate (zqdt(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,             &
                                      Mgrid%kbd:Mgrid%ked))

  tdt_phys = 0. ; qdt_phys = 0. ; sdt_phys = 0. ; sdt_vdif = 0. ; sdt_drag = 0.
  tdt_mix = 0. ; qdt_mix = 0.
  udt_mix = 0. ; vdt_mix = 0.

  tdt_damp = 0.
  qdt_damp = 0.
  udt_damp = 0.
  vdt_damp = 0.

!-----------------------------------------------------------------------
! allocate smooth arrays
!-----------------------------------------------------------------------

  allocate (umean(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,             &
                                       Mgrid%kbd:Mgrid%ked))
  allocate (vmean(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,             &
                                       Mgrid%kbd:Mgrid%ked))
  allocate (tmean(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,             &
                                       Mgrid%kbd:Mgrid%ked))
  allocate (qmean(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,             &
                                       Mgrid%kbd:Mgrid%ked))
  allocate (udiv(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,             &
                                      Mgrid%kbd:Mgrid%ked))
  allocate (vdiv(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,             &
                                      Mgrid%kbd:Mgrid%ked))
  allocate (var(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                                     Mgrid%kbd:Mgrid%ked))

!-----------------------------------------------------------------------
! allocate tracer tendency
!-----------------------------------------------------------------------

  allocate (f_qc(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,             &
                                      Mgrid%kbd:Mgrid%ked))
  allocate (uatm(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,             &
                                      Mgrid%kbd:Mgrid%ked))
  allocate (vatm(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,             &
                                      Mgrid%kbd:Mgrid%ked))

!-----------------------------------------------------------------------
! allocate diagnostic arrays
!-----------------------------------------------------------------------

  allocate (rain(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed))
  allocate (snow(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed))

  rain = 0. ; snow = 0.

!-----------------------------------------------------------------------
! allocate diffusion arrays
!-----------------------------------------------------------------------

  allocate (diffu(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,             &
                                       Mgrid%kbd:Mgrid%ked))
  allocate (diffv(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,             &
                                       Mgrid%kbd:Mgrid%ked))
  allocate (diffh(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,             &
                                       Mgrid%kbd:Mgrid%ked))
  allocate (cdu  (Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed))
  allocate (cdv  (Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed))
  allocate (cdh  (Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed))
  allocate (dflxu(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed))
  allocate (dflxv(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed))
  allocate (dflxh(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed))
  allocate (vsatu(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed))
  allocate (vsatv(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed))
  allocate (vsatm(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed))

!-----------------------------------------------------------------------
! initialize time pointers and get timestep and runtime
!-----------------------------------------------------------------------

  call get_time_pointers ( past, pres, future )

  call time_manager_init

  call get_time ( Time_last - Time_init, secs )
  run_time = real(secs)

!-----------------------------------------------------------------------
! initialize ncdf_io
!-----------------------------------------------------------------------

  call ncdf_io_init ( Mgrid, Time_step )

!-----------------------------------------------------------------------
! initialize restart read
!-----------------------------------------------------------------------

  call read_restart_init (                                             &
                                               Mgrid, Time, Time_step, &
                                   rlonmin, rlonmax, rlatmin, rlatmax, &
                                                                zetaw, &
                                                         do_coldstart )

!-----------------------------------------------------------------------
! read restart
!-----------------------------------------------------------------------

  call read_restart (                                                  &
                                                                Mgrid, &
                                                           topog, sst, &
                                           ppref, uuref, taref, qvref, &
                                      ps, dpm, uu, vv, oo, ta, qv, qc, &
                 udt_filt, vdt_filt, ttend, qtend, tdt_filt, qdt_filt, &
                               umean, vmean, tmean, qmean, udiv, vdiv, &
                                                           past, pres )

!-----------------------------------------------------------------------
! define horizontal metric terms and horizontal metric derived-type 
!-----------------------------------------------------------------------

  call horiz_metric (                                                  &
                                                       Mgrid, Hmetric, &
                                   rlonmin, rlonmax, rlatmin, rlatmax, &
                                         do_cartesian, do_cylindrical )

!-----------------------------------------------------------------------
! read topography
!-----------------------------------------------------------------------

!stg  call set_topog ( Mgrid, Hmetric, topog, sst )

!-----------------------------------------------------------------------
! define vertical metric terms and vertical metric derived-type
!-----------------------------------------------------------------------

  call vert_metric (                                                   &
                                              Mgrid, Hmetric, Vmetric, &
                                                    pbot, ptop, zetaw, &
                                                                topog )

!-----------------------------------------------------------------------
! initializer restart write
!-----------------------------------------------------------------------

  call write_restart_init ( Mgrid, Hmetric, Vmetric )

!-----------------------------------------------------------------------
! initialize history writes
!-----------------------------------------------------------------------

  call write_history_init ( Mgrid, Hmetric, Vmetric, Time, Time_step )

!-----------------------------------------------------------------------
! register diagnostics
!-----------------------------------------------------------------------

  axid_mm = (/axid(xm_),axid(ym_)/)
  axid_mmm = (/axid(xm_),axid(ym_),axid(zm_)/)

  name = "tdt_vdif"
  longname = "temperature tendency from vertical diffusion"
  units = 'K/s'
  id_tdt_vdif = register_diag_field (                                  &
                      model, name, axid_mmm, Time,                     &
                      longname, units, silly, range=(/-50., 50./) )

  name = "qdt_vdif"
  longname = "specific humidity tendency from vertical diffusion"
  units = 'kg/kg/s'
  id_qdt_vdif = register_diag_field (                                  &
                      model, name, axid_mmm, Time,                     &
                      longname, units, silly, range=(/-50., 50./) )
  name = "tdt_phys"
  longname = "temperature tendency from physics"
  units = 'K/s'
  id_tdt_phys = register_diag_field (                                  &
                      model, name, axid_mmm, Time,                     &
                      longname, units, silly, range=(/-10., 10./) )
  name = "qdt_phys"
  longname = "specific humidity tendency from physics"
  units = '1/s'
  id_qdt_phys = register_diag_field (                                  &
                      model, name, axid_mmm, Time,                     &
                      longname, units, silly, range=(/-0.1, 0.1/) )
  name = "tdt_damp"
  longname = "temperature tendency from damping"
  units = 'K/s'
  id_tdt_damp = register_diag_field (                                  &
                      model, name, axid_mmm, Time,                     &
                      longname, units, silly, range=(/-0.1, 0.1/) )
  name = "qdt_damp"
  longname = "specific humidity tendency from damping"
  units = '1/s'
  id_qdt_damp = register_diag_field (                                  &
                      model, name, axid_mmm, Time,                     &
                      longname, units, silly, range=(/-0.02, 0.02/) )
  name = "sdt_phys"
  longname = "entropy tendency from physics"
  units = '1/s'
  id_sdt_phys = register_diag_field (                                  &
                      model, name, axid_mmm, Time,                     &
                      longname, units, silly, range=(/-0.1, 0.1/) )
  name = "sdt_vdif"
  longname = "entropy tendency from diffusion"
  units = '1/s'
  id_sdt_phys = register_diag_field (                                  &
                      model, name, axid_mmm, Time,                     &
                      longname, units, silly, range=(/-0.1, 0.1/) )
  name = "sdt_diss"
  longname = "entropy tendency from dissipation"
  units = 'm2/s3/K'
  id_sdt_diss = register_diag_field (                                  &
                      model, name, axid_mmm, Time,                     &
                      longname, units, silly, range=(/-0.1, 0.1/) )
  name = "sdt_adv"
  longname = "entropy tendency from advection errors"
  units = 'm2/s3/K'
  id_sdt_adv = register_diag_field (                                   &
                      model, name, axid_mmm, Time,                     &
                      longname, units, silly, range=(/-0.3, 0.3/) )
  name = "sdt_dry"
  longname = "entropy tendency from evaporation"
  units = 'm2/s3/K'
  id_sdt_dry = register_diag_field (                                   &
                      model, name, axid_mmm, Time,                     &
                      longname, units, silly, range=(/-0.1, 0.1/) )
  name = "sdt_drag"
  longname = "entropy tendency from condensate drag"
  units = 'm2/s3/K'
  id_sdt_drag = register_diag_field (                                  &
                      model, name, axid_mmm, Time,                     &
                      longname, units, silly, range=(/-0.1, 0.1/) )
  name = "tdt_cond"
  longname = "microphysics temperature tendency"
  units = 'K/s'
  id_tdt_cond = register_diag_field (                                  &
                      model, name, axid_mmm, Time,                     &
                      longname, units, silly, range=(/-10., 10./) )
  name = "qdt_cond"
  longname = "microphysics specific humidity tendency"
  units = '1/s'
  id_qdt_cond = register_diag_field (                                  &
                      model, name, axid_mmm, Time,                     &
                      longname, units, silly, range=(/-0.1, 0.1/) )
  name = "qdt_drag"
  longname = "energy tendency from condensate drag"
  units = 'm2/s3'
  id_qdt_drag = register_diag_field (                                  &
                       model, name, axid_mmm, Time,                    &
                       longname, units, silly, range=(/-1.e4, 1.e4/) )
  name = "alfomega"
  longname = "alfa*dp/dt"
  units = 'm2/s3'
  id_alfomega = register_diag_field (                                  &
                       model, name, axid_mmm, Time,                    &
                       longname, units, silly, range=(/-1.e4, 1.e4/) )
  name = "tdt_vdif_int"
  longname = "column sensible heating from vertical diffusion"
  units = 'W/m2'
  id_tdt_vdif_int = register_diag_field (                              &
                       model, name, axid_mm, Time,                     &
                       longname, units, silly, range=(/-1.e7, 1.e7/) )
  name = "qdt_vdif_int"
  longname = "column latent heating from vertical diffusion"
  units = 'W/m2'
  id_qdt_vdif_int = register_diag_field (                              &
                       model, name, axid_mm, Time,                     &
                       longname, units, silly, range=(/-3.e6, 3.e6/) )
  name = "tdt_phys_int"
  longname = "column heating from physics"
  units = 'W/m2'
  id_tdt_phys_int = register_diag_field (                              &
                       model, name, axid_mm, Time,                     &
                       longname, units, silly, range=(/-2.e6, 2.e6/) )
  name = "qdt_phys_int"
  longname = "column moistening from physics"
  units = 'W/m2'
  id_qdt_phys_int = register_diag_field (                              &
                       model, name, axid_mm, Time,                     &
                       longname, units, silly, range=(/-1.e9, 1.e9/) )
  name = "tdt_damp_int"
  longname = "column heating from damping"
  units = 'W/m2'
  id_tdt_damp_int = register_diag_field (                              &
                       model, name, axid_mm, Time,                     &
                       longname, units, silly, range=(/-1.e4, 1.e4/) )
  name = "qdt_damp_int"
  longname = "column moistening from damping"
  units = 'W/m2'
  id_qdt_damp_int = register_diag_field (                              &
                       model, name, axid_mm, Time,                     &
                       longname, units, silly, range=(/-1.e9, 1.e9/) )
  name = "qdt_drag_int"
  longname = "column entropy tendency from condensate drag"
  units = 'W/m2'
  id_qdt_drag_int = register_diag_field (                              &
                       model, name, axid_mm, Time,                     &
                       longname, units, silly, range=(/-1.e4, 1.e4/) )
  name = "sdt_phys_int"
  longname = "column entropy tendency from physics"
  units = 'W/m2/K'
  id_sdt_phys_int = register_diag_field (                              &
                       model, name, axid_mm, Time,                     &
                       longname, units, silly, range=(/-2.e5, 2.e5/) )
  name = "sdt_vdif_int"
  longname = "column entropy tendency from diffusion"
  units = 'W/m2/K'
  id_sdt_vdif_int = register_diag_field (                              &
                       model, name, axid_mm, Time,                     &
                       longname, units, silly, range=(/-2.e5, 2.e5/) )
  name = "sdt_diss_int"
  longname = "column entropy tendency from dissipation"
  units = 'W/m2/K'
  id_sdt_diss_int = register_diag_field (                              &
                       model, name, axid_mm, Time,                     &
                       longname, units, silly, range=(/-1.e4, 1.e4/) )
  name = "sdt_adv_int"
  longname = "column entropy tendency from advection errors"
  units = 'W/m2/K'
  id_sdt_adv_int = register_diag_field (                               &
                       model, name, axid_mm, Time,                     &
                       longname, units, silly, range=(/-1.e5, 1.e5/) )
  name = "sdt_dry_int"
  longname = "column entropy tendency from evaporation"
  units = 'W/m2/K'
  id_sdt_dry_int = register_diag_field (                               &
                       model, name, axid_mm, Time,                     &
                       longname, units, silly, range=(/-1.e5, 1.e5/) )
  name = "sdt_drag_int"
  longname = "column entropy tendency from condensate drag"
  units = 'W/m2/K'
  id_sdt_drag_int = register_diag_field (                              &
                       model, name, axid_mm, Time,                     &
                       longname, units, silly, range=(/-1.e4, 1.e4/) )
  name = "tdt_cond_int"
  longname = "column heating from microphsysics"
  units = 'W/m2'
  id_tdt_cond_int = register_diag_field (                              &
                       model, name, axid_mm, Time,                     &
                       longname, units, silly, range=(/-1.e6, 1.e6/) )
  name = "qdt_cond_int"
  longname = "column moistening from microphysics"
  units = 'W/m2'
  id_qdt_cond_int = register_diag_field (                              &
                       model, name, axid_mm, Time,                     &
                       longname, units, silly, range=(/-1.e4, 1.e4/) )

!-----------------------------------------------------------------------
! get base time specified in diag table
!-----------------------------------------------------------------------

  call get_base_date ( date(1), date(2), date(3),                      &
                       date(4), date(5), date(6) )

  Time_base = set_date ( date(1), date(2), date(3),                    &
                         date(4), date(5), date(6) )

!-----------------------------------------------------------------------
! initialize core
!-----------------------------------------------------------------------

  call zetac_init (                                                    &
                                              Mgrid, Hmetric, Vmetric, &
                                           Time, Time_step, Time_base, &
                                    ppref, uuref, taref, qvref, gzref, &
                                        qc_index, nsteps, ndt_forward )

!-----------------------------------------------------------------------
! initialize large-scale thermal forcing
!-----------------------------------------------------------------------

  call ls_forcing_init ( Mgrid, Vmetric, Time )

!-----------------------------------------------------------------------
! initialize microphysics ('num_trac' will be incremented)
!-----------------------------------------------------------------------

  call moisture_driver_init ( Mgrid, Time, qvmin )

!-----------------------------------------------------------------------
! is the first timestep a forward step?
!-----------------------------------------------------------------------

  if ( ndt_forward > 0 ) then
     call get_time ( Time_step, secs )
     secs_forward = ndt_forward*secs
     call get_time ( Time - Time_base, secs )
     do_forward = ( mod(secs,secs_forward) == 0 )
  endif

  call mpp_clock_end ( atmosphere_init_clock )

  return
end subroutine atmosphere_init

!#######################################################################

subroutine atmosphere_end ( Time  )

type (time_type),                                   intent(in) ::      &
                                                                  Time

real :: time_now, time_end

  call mpp_clock_begin ( atmosphere_end_clock )
  
!-----------------------------------------------------------------------
! finalize core
!-----------------------------------------------------------------------

  call zetac_end
 
!-----------------------------------------------------------------------
! reset domain decomposition for output
!-----------------------------------------------------------------------

  call set_domain ( Mgrid%Domain_nohalo )

!-----------------------------------------------------------------------
! final times
!-----------------------------------------------------------------------

  call get_time (Time - Time_step - Time_base, secs)
  time_now = real(secs)/60.0          ! restart time axis is in minutes
  call get_time (Time - Time_base, secs)
  time_end = real(secs)/60.0

  call get_time_pointers ( past, pres, future )

!-----------------------------------------------------------------------
! finalize moisture
!-----------------------------------------------------------------------

  call moisture_driver_end ( Mgrid, time_now, past )
  call moisture_driver_end ( Mgrid, time_end, pres )

!-----------------------------------------------------------------------
! write restart
!-----------------------------------------------------------------------

  call update_halos ( Mgrid, ttend )
  call update_halos ( Mgrid, qtend )

  call write_restart (                                                 &
                                               Mgrid, Time, Time_step, &
                                                           topog, sst, &
                                           ppref, uuref, taref, qvref, &
                                      ps, dpm, uu, vv, oo, ta, qv, qc, &
                 udt_filt, vdt_filt, ttend, qtend, tdt_filt, qdt_filt, &
                               umean, vmean, tmean, qmean, udiv, vdiv, &
                                                           past, pres )
					
  call mpp_clock_end ( atmosphere_end_clock )

  return
end subroutine atmosphere_end

!#######################################################################

end module atmosphere_mod
