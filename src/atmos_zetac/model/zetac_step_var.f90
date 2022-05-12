module zetac_step_var_mod

!-----------------------------------------------------------------------
! module information
!-----------------------------------------------------------------------

use mpp_mod,            only : mpp_pe, mpp_root_pe
use fms_mod,            only : write_version_number,                   &
                               file_exist, close_file, stdlog, stdout
use diag_manager_mod,   only : register_diag_field
use field_manager_mod,  only : model_atmos
use time_manager_mod,   only : time_type
use tracer_manager_mod, only : get_tracer_names, get_tracer_index
use vert_advection_mod, only : FLUX_FORM

use zetac_axis_names_mod,        only : xm_, ym_, zm_, xu_, yu_,       &
                                        xv_, yv_, zw_ 
use zetac_axes_mod,              only : axid,                          &
                                        gridu, gridv, gridm, gridw
use zetac_cgrid_interp_mod,      only : cgrid_interp_uv
use zetac_extrap_var_mod,        only : extrap_ptq, extrap_uv,         &
                                        extrap_tracer
use zetac_horiz_grid_type_mod,   only : horiz_grid_type
use zetac_horiz_metric_type_mod, only : horiz_metric_type
use zetac_ncdf_io_mod,           only : ncdf_fms_write
use zetac_tracer_mod,            only : get_method
use zetac_time_pointers_mod,     only : ntime
use zetac_tridiag_mod,           only : solve_tridiag
use zetac_update_halos_mod,      only : update_halos
use zetac_vert_diff_mod,         only : vert_diff
use zetac_vert_metric_type_mod,  only : vert_metric_type

!-----------------------------------------------------------------------
! public interfaces
!-----------------------------------------------------------------------

implicit none
private
public step_ps, step_uv, step_tq, step_tracer, step_var_init

real, allocatable, dimension(:,:,:) :: g_uu, g_vv, g_ta, g_qv
real, allocatable, dimension(:,:)   :: uref, tref, qref
real, allocatable, dimension(:)     :: pref
real, allocatable, dimension(:)     :: dxu, cosulat, cosvlat

character (len=*),  parameter :: module='zetac_step_var_mod'
character (len=16), parameter :: model='atmos_mod'
character (len=16) :: name, units
character (len=64) :: longname

real, parameter :: silly=-9999.9

integer :: ibc, iec, ibd, ied
integer :: jbc, jec, jbd, jed
integer :: kbd, ked

!-----------------------------------------------------------------------
! version control information 
!-----------------------------------------------------------------------

character(len=128)  :: version =  '$Id: zetac_step_var.f90,v 1.1.2.9.2.15 2005/07/30 02:45:23 stg Exp $'
character(len=128)  :: tag     =  '$Name:  $'

contains

!#######################################################################

subroutine step_ps (                                                   &
                                                                Mgrid, &
                                                             ps, f_ps, &
                                                         past, future, &
                                                                 delt )

type (horiz_grid_type),                     intent (in)    ::          &
                                                                Mgrid

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed, ntime),      &
                                            intent (inout) ::          &
                                                                   ps

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed),             &
                                            intent (in)    ::          &
                                                                 f_ps

real,                                       intent (in)    ::          &
                                                                 delt

integer, intent (in) :: past, future

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

integer :: i, j

!-----------------------------------------------------------------------
! explicit update of surface pressure
!-----------------------------------------------------------------------

  do j=jbc,jec
     do i=ibc,iec
        ps(i,j,future) = ps(i,j,past) + delt*f_ps(i,j)
     enddo
  enddo

  call update_halos ( Mgrid, ps(:,:,future) )

  return
end subroutine step_ps

!#######################################################################

subroutine step_tq (                                                   &
                                                     Mgrid, Time_next, &
                                                       ps, ta, qv, dp, &
                                                           f_ta, f_qv, &
                                                   tdt_vdif, qdt_vdif, &
                                                     beta, gam, ldiag, &
                                       do_flux, past, present, future, &
                                                                 delt )

type (horiz_grid_type),                      intent (in)    ::         &
                                                                Mgrid

type (time_type),                            intent (in)    ::         &
                                                            Time_next

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed, ntime),      &
                                             intent (inout) ::         &
                                                                   ps

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked, ntime), intent (inout) ::         &
                                                               ta, qv

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked, ntime), intent (in)    ::         &
                                                                   dp

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),        intent (in)    ::         &
                                                           f_ta, f_qv

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),        intent (inout) ::         &
                                                   tdt_vdif, qdt_vdif

real, dimension(Mgrid%ibc:Mgrid%iec, Mgrid%jbc:Mgrid%jec,              &
                Mgrid%kbd+1:Mgrid%ked),      intent (inout) ::         &
                                                     beta, gam, ldiag

real,                                        intent (in)    ::         &
                                                                 delt

logical, intent (in) :: do_flux
integer, intent (in) :: past, present, future

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

integer :: i, j, k

!-----------------------------------------------------------------------
! explicit updates of ta and qv
!-----------------------------------------------------------------------

  if (do_flux) then

     do k=kbd+1,ked
        do j=jbc,jec
           do i=ibc,iec
              g_ta(i,j,k) = (dp(i,j,k,past)*ta(i,j,k,past)             &
                 + delt*dp(i,j,k,present)*f_ta(i,j,k))/dp(i,j,k,future)
              g_qv(i,j,k) = (dp(i,j,k,past)*qv(i,j,k,past)             &
                 + delt*dp(i,j,k,present)*f_qv(i,j,k))/dp(i,j,k,future)
           enddo
        enddo
     enddo

  else

     do k=kbd+1,ked
        do j=jbc,jec
           do i=ibc,iec
              g_ta(i,j,k) = ta(i,j,k,past) + delt*f_ta(i,j,k)
              g_qv(i,j,k) = qv(i,j,k,past) + delt*f_qv(i,j,k)
           enddo
        enddo
     enddo

  endif

!-----------------------------------------------------------------------
! solve tridiagonal system for ta and qv
!-----------------------------------------------------------------------

  call vert_diff (                                                     &
                                                                Mgrid, &
                                               ta(:,:,:,future), g_ta, &
                                                     beta, gam, ldiag )
  call vert_diff (                                                     &
                                                                Mgrid, &
                                               qv(:,:,:,future), g_qv, &
                                                     beta, gam, ldiag )

!-----------------------------------------------------------------------
! extrapolate outside open boundaries
!-----------------------------------------------------------------------

  call extrap_ptq (                                                    &
                                                                Mgrid, &
                                                           ps, ta, qv, &
                                                     pref, tref, qref, &
                                                                 delt )

!-----------------------------------------------------------------------
! diagnose total damping and diffusion
!-----------------------------------------------------------------------

  if (do_flux) then

     do k=kbd+1,ked
        tdt_vdif(:,:,k) = tdt_vdif(:,:,k)                              &
          + (ta(:,:,k,future) - g_ta(:,:,k))*dp(:,:,k,future)/         &
                                       (delt*dp(:,:,k,present))
        qdt_vdif(:,:,k) = qdt_vdif(:,:,k)                              &
          + (qv(:,:,k,future) - g_qv(:,:,k))*dp(:,:,k,future)/         &
                                       (delt*dp(:,:,k,present))
     enddo

  else

     do k=kbd+1,ked
        tdt_vdif(:,:,k) = tdt_vdif(:,:,k)                              &
                             + (ta(:,:,k,future) - g_ta(:,:,k))/delt
        qdt_vdif(:,:,k) = qdt_vdif(:,:,k)                              &
                             + (qv(:,:,k,future) - g_qv(:,:,k))/delt
     enddo

  endif

  tdt_vdif(:,:,kbd) = 0.
  qdt_vdif(:,:,kbd) = 0.

!-----------------------------------------------------------------------
! update halos for next fast timestep
!-----------------------------------------------------------------------

  call update_halos ( Mgrid, ta(:,:,:,future) )
  call update_halos ( Mgrid, qv(:,:,:,future) )

  return
end subroutine step_tq

!#######################################################################

subroutine step_uv (                                                   &
                                                     Mgrid, Time_next, &
                                                     uu, vv, dpu, dpv, &
                                                           f_uu, f_vv, &
                                                   udt_vdif, vdt_vdif, &
                                             beta_u, gamma_u, ldiag_u, &
                                             beta_v, gamma_v, ldiag_v, &
                                       do_flux, past, present, future, &
                                                                 delt )

type (horiz_grid_type),                     intent (in)    ::          &
                                                                Mgrid

type (time_type),                           intent (in)    ::          &
                                                            Time_next

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked, ntime), intent (inout) ::         &
                                                               uu, vv

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked, ntime), intent (in)    ::         &
                                                             dpu, dpv

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),        intent (in)    ::         &
                                                           f_uu, f_vv

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),       intent (inout) ::          &
                                                   udt_vdif, vdt_vdif

real, dimension(Mgrid%ibc:Mgrid%iec, Mgrid%jbc:Mgrid%jec,              &
                Mgrid%kbd+1:Mgrid%ked),     intent (inout) ::          &
                   beta_u, gamma_u, ldiag_u, beta_v, gamma_v, ldiag_v

real,                                       intent (in)    ::          &
                                                                 delt

logical, intent (in) :: do_flux
integer, intent (in) :: past, present, future

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

integer :: i, j, k

!-----------------------------------------------------------------------
! explicit updates of uu and vv
!-----------------------------------------------------------------------

  if (do_flux) then

     do k=kbd+1,ked
        do j=jbc,jec
           do i=ibc,iec
              g_uu(i,j,k) = (dpu(i,j,k,past)*uu(i,j,k,past)            &
               + delt*dpu(i,j,k,present)*f_uu(i,j,k))/dpu(i,j,k,future)
              g_vv(i,j,k) = (dpv(i,j,k,past)*vv(i,j,k,past)            &
               + delt*dpv(i,j,k,present)*f_vv(i,j,k))/dpv(i,j,k,future)
           enddo
        enddo
     enddo

  else

     do k=kbd+1,ked
        do j=jbc,jec
           do i=ibc,iec
              g_uu(i,j,k) = uu(i,j,k,past) + delt*f_uu(i,j,k)
              g_vv(i,j,k) = vv(i,j,k,past) + delt*f_vv(i,j,k)
           enddo
        enddo
     enddo

  endif

!-----------------------------------------------------------------------
! solve tridiagonal system for uu and vv
!-----------------------------------------------------------------------

  call vert_diff (                                                     &
                                                                Mgrid, &
                                               uu(:,:,:,future), g_uu, &
                                             beta_u, gamma_u, ldiag_u )
  call vert_diff (                                                     &
                                                                Mgrid, &
                                               vv(:,:,:,future), g_vv, &
                                             beta_v, gamma_v, ldiag_v )

  call extrap_uv ( Mgrid, uu(:,:,:,future), vv(:,:,:,future), uref )

!-----------------------------------------------------------------------
! diagnose total damping and diffusion
!-----------------------------------------------------------------------

  if (do_flux) then

     do k=kbd+1,ked
        udt_vdif(:,:,k) = udt_vdif(:,:,k)                              &
              + (uu(:,:,k,future) - g_uu(:,:,k))*dpu(:,:,k,future)/    &
                                           (delt*dpu(:,:,k,present))
        vdt_vdif(:,:,k) = vdt_vdif(:,:,k)                              &
              + (vv(:,:,k,future) - g_vv(:,:,k))*dpv(:,:,k,future)/    &
                                           (delt*dpv(:,:,k,present))
     enddo

  else

     do k=kbd+1,ked
        udt_vdif(:,:,k) = udt_vdif(:,:,k)                              &
                             + (uu(:,:,k,future) - g_uu(:,:,k))/delt
        vdt_vdif(:,:,k) = vdt_vdif(:,:,k)                              &
                             + (vv(:,:,k,future) - g_vv(:,:,k))/delt
     enddo

  endif

  udt_vdif(:,:,kbd) = 0.
  vdt_vdif(:,:,kbd) = 0.

!-----------------------------------------------------------------------
! update halos
!-----------------------------------------------------------------------

  call update_halos ( Mgrid, uu(:,:,:,future) )
  call update_halos ( Mgrid, vv(:,:,:,future) )

  return
end subroutine step_uv

!#######################################################################

subroutine step_tracer (                                               &
                                                                Mgrid, &
                                                             tr, f_tr, &
                                                         past, future, &
                                                                 delt )

type (horiz_grid_type),                     intent (in)    ::          &
                                                                Mgrid

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked, ntime),                           &
                                            intent (inout) ::          &
                                                                   tr

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),                                  &
                                            intent (in)    ::          &
                                                                 f_tr

real,                                       intent (in)    ::          &
                                                                 delt

integer,                                    intent (in)    ::          &
                                                         past, future

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

integer :: num_form
integer :: i, j, k

!-----------------------------------------------------------------------
! update tracer
!-----------------------------------------------------------------------

  do k=kbd+1,ked
     do j=jbc,jec
        do i=ibc,iec
           tr(i,j,k,future) = tr(i,j,k,past) + delt*f_tr(i,j,k)
        enddo
     enddo
  enddo

!  call extrap_tracer ( Mgrid, tr(:,:,:,future), delt )

!-----------------------------------------------------------------------
! update halos
!-----------------------------------------------------------------------

!  call update_halos ( Mgrid, tr(:,:,:,future) )

  return
end subroutine step_tracer

!#######################################################################

subroutine step_var_init (                                             &
                                        Mgrid, Hmetric, Vmetric, Time, &
                                           ppref, uuref, taref, qvref )

type (horiz_grid_type),                         intent (in)   ::       &
                                                                Mgrid

type (horiz_metric_type),                       intent (in)   ::       &
                                                              Hmetric

type (vert_metric_type),                        intent (in)   ::       &
                                                              Vmetric

type (time_type),                               intent (in)   ::       &
                                                                 Time

real, dimension(Mgrid%jbg:Mgrid%jeg, Mgrid%iz),                        &
                                            intent (in)    ::          &
                                           ppref, uuref, taref, qvref

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

integer :: io, ierr
integer :: i, j, k
 
!-----------------------------------------------------------------------
! write version number to logfile
!-----------------------------------------------------------------------

  call write_version_number (version, tag)

  ibc = Mgrid%ibc
  iec = Mgrid%iec
  ibd = Mgrid%ibd
  ied = Mgrid%ied

  jbc = Mgrid%jbc
  jec = Mgrid%jec
  jbd = Mgrid%jbd
  jed = Mgrid%jed

  kbd = Mgrid%kbd
  ked = Mgrid%ked

!-----------------------------------------------------------------------
! arrays for tridiagonal problem
!-----------------------------------------------------------------------

  allocate (g_uu(ibd:ied,jbd:jed,kbd:ked))
  allocate (g_vv(ibd:ied,jbd:jed,kbd:ked))
  allocate (g_ta(ibd:ied,jbd:jed,kbd:ked))
  allocate (g_qv(ibd:ied,jbd:jed,kbd:ked))

  g_uu = 0.
  g_vv = 0.
  g_ta = 0.
  g_qv = 0.

  allocate (pref(jbd:jed))
  allocate (uref(jbd:jed,kbd:ked))
  allocate (tref(jbd:jed,kbd:ked))
  allocate (qref(jbd:jed,kbd:ked))

  pref = ppref(jbd:jed,ked)
  uref = uuref(jbd:jed,:)
  tref = taref(jbd:jed,:)
  qref = qvref(jbd:jed,:)

  return
end subroutine step_var_init

!#######################################################################

end module zetac_step_var_mod

