module zetac_vert_turb_mod

!-----------------------------------------------------------------------
! module information1
!-----------------------------------------------------------------------

use fms_mod,          only : write_version_number, mpp_pe, mpp_root_pe,&
                             check_nml_error, stdlog, error_mesg, FATAL
use mpp_mod,          only : input_nml_file
use mpp_domains_mod,  only : mpp_update_domains, eupdate, nupdate

use constants_mod,    only : Rdgas, Grav

use diag_manager_mod,            only : register_diag_field
use field_manager_mod,           only : model_atmos
use time_manager_mod,            only : time_type

use zetac_axes_mod,              only : axid, gridm
use zetac_axis_names_mod,        only : xm_, ym_
use zetac_cgrid_interp_mod,      only : cgrid_interp_uv,               &
                                        cgrid_interp_mass
use zetac_phys_con_mod,          only : rrat
use zetac_horiz_grid_type_mod,   only : horiz_grid_type
use zetac_horiz_metric_type_mod, only : horiz_metric_type
use zetac_ncdf_io_mod,           only : ncdf_fms_write
use zetac_time_pointers_mod,     only : ntime
use zetac_tridiag_mod,           only : solve_tridiag
use zetac_update_halos_mod,      only : update_halos
use zetac_vert_metric_type_mod,  only : vert_metric_type
use zetac_update_halos_mod,      only : update_halos

!-----------------------------------------------------------------------
! public interfaces
!-----------------------------------------------------------------------

implicit none
private
public vert_turb, vert_turb_init

character (len=*),  parameter :: module='zetac_vert_turb_mod'
character (len=16), parameter :: model='atmos_mod'
character (len=16) :: name, units
character (len=64) :: longname

real :: ptop

! namelist:
real :: diffm_const=0., diffh_const=0.
real :: cdm_const=0.001, cdh_const=0.001
real :: z0=1.0, zh=1.0
real :: Tlapse, qlapse
real :: gust=1.0
real :: shap=0.0
real :: tau_diff=0.0
logical :: conserve_etot, conserve_mass, conserve_eh2o

namelist /vert_turb_nml/         &  !note: this namelist also read
  diffm_const, diffh_const,      &  !      by calc_var_tend
  cdm_const, cdh_const,          &
  z0, zh,                        &
  Tlapse, qlapse, gust, shap,    &
  tau_diff,                      &
  conserve_etot,                 &
  conserve_mass,                 &
  conserve_eh2o

real, allocatable, dimension(:)   :: logz1z0, logz1zh, pp, vs, zz
real, allocatable, dimension(:)   :: zetam, zetaw, dzetam, dzetaw
real, allocatable, dimension(:,:) :: dzref
real, allocatable, dimension(:,:) :: uatm, vatm, v2, v2atu, v2atv, cdm
real, allocatable, dimension(:,:) :: v2old

real, parameter :: tiny = 1.e-9
real, parameter :: silly=-9999.9

real :: tau = -1.
real :: fac

integer :: i, j, k
integer :: ibc, iec, ibd, ied
integer :: jbc, jec, jbd, jed
integer :: kbd, ked, iz

integer :: id_cdm, id_cdh

!-----------------------------------------------------------------------
! version control information
!-----------------------------------------------------------------------

character(len=128)  :: version =  '$Id: zetac_vert_turb.f90,v 1.1.2.9.2.4.2.5 2005/07/18 18:44:25 stg Exp $'
character(len=128)  :: tag     =  '$Name:  $'

contains

!#######################################################################

subroutine vert_turb (                                                 &
                                                          Mgrid, Time, &
                                                               uu, vv, &
                                                           ps, ta, qv, &
                                                  diffu, diffv, diffh, &
              cdu, cdv, cdh, dflxu, dflxv, dflxh, vsatu, vsatv, vsatm, &
                                                                 delt )

type (horiz_grid_type),                      intent (inout) ::         &
                                                                Mgrid
       
type (time_type),                            intent (in)    ::         &
                                                                 Time

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed),             &
                                             intent (in)    ::         &
                                                               uu, vv

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed),             &
                                             intent (in)    ::         &
                                                           ps, ta, qv

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),        intent (out)   ::         &
                                                  diffu, diffv, diffh

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed),             &
                                             intent (out)   ::         &
                                                        cdu, cdv, cdh

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed),             &
                                             intent (out)   ::         &
                             dflxu, dflxv, dflxh, vsatu, vsatv, vsatm

real, intent(in) :: delt

  diffu = 0.
  diffv = 0.
  diffh = 0.
  do k=kbd+1,ked-1
     if (zetaw(k) > 0.10) then
        fac = (dzetaw(k)/dzetaw(iz-1))**2
        do j=jbd,jed
           diffh(:,j,k) = diffh_const * fac
        enddo
     endif
  enddo
  do k=kbd+1,ked-1
     if (zetaw(k) > 0.10) then
        fac = (dzetaw(k)/dzetaw(iz-1))**2
        do j=jbd,jed
           diffu(:,j,k) = diffm_const * fac
           diffv(:,j,k) = diffm_const * fac
        enddo
     endif
  enddo

  do k=kbd+1,ked-1
     fac = (zetaw(k)/zetaw(ked-1))!**2
     diffh(:,:,k) = diffh(:,:,k)*fac
     diffu(:,:,k) = diffu(:,:,k)*fac
     diffv(:,:,k) = diffv(:,:,k)*fac
  enddo

  call cgrid_interp_mass ( Mgrid, uu, vv, uatm, vatm )

  v2 = uatm*uatm + vatm*vatm + gust*gust

  if (tau >= 0.) then
     tau = dim(1., delt/(tau_diff + tiny))
  else
     tau = 0.
  endif
  v2 = tau*v2old + (1. - tau)*v2
  v2old = v2
  vsatm = sqrt(v2)

  do i=ibd,ied
     pp(:) = ptop + (ps(i,:) - ptop)*zetam(ked)
     zz(:) = Rdgas*ta(i,:)*(1.0 + qv(i,:)*rrat)/Grav*log(ps(i,:)/pp(:))
     vs(:) = (logz1z0(:) / log(zz(:)/z0))**2
     cdm(i,:) = (cdm_const + 0.001*shap*vsatm(i,:))* vs(:)
     vs(:) = (logz1zh(:)*logz1z0(:)) / (log(zz(:)/zh)*log(zz(:)/z0))
     cdh(i,:) = cdh_const* vs(:)
  enddo

  call cgrid_interp_uv ( Mgrid, cdm, cdu, cdv )

!!$  dflxh = vsatm          ! for implicit vertical diffusion
  dflxh = 0.

  call mpp_update_domains ( v2, Mgrid%Domain, flags=eupdate+nupdate )
  call cgrid_interp_uv ( Mgrid, v2, v2atu, v2atv )
  vsatu = sqrt(v2atu)
  vsatv = sqrt(v2atv)
!!$  dflxu = (v2atu + uu*uu + gust*gust)/(vsatu + tiny)
!!$  dflxv = (v2atv + vv*vv + gust*gust)/(vsatv + tiny)
dflxu=0.
dflxv=0.

  call ncdf_fms_write (Mgrid, Time, id_cdm, gridm, cdm)
  call ncdf_fms_write (Mgrid, Time, id_cdh, gridm, cdh)

  return
end subroutine vert_turb

!#######################################################################

subroutine vert_turb_init ( Mgrid, Vmetric, Time, zmref, zref, pref )

type (horiz_grid_type),                     intent (in)  ::            &
                                                                 Mgrid

type (vert_metric_type),                     intent (in)    ::         &
                                                               Vmetric

type (time_type),                           intent (in)  ::            &
                                                                  Time

real, dimension(Mgrid%jbg:Mgrid%jeg),       intent (in)  ::      zmref

real, dimension(Mgrid%jbg:Mgrid%jeg,Mgrid%iz),                         &
                                            intent (in)  ::       zref

real, dimension(Mgrid%jbg:Mgrid%jeg,Mgrid%iz),                         &
                                            intent (in)  ::       pref

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

integer, dimension(2) :: axid_mm
integer :: io, ierr

!-----------------------------------------------------------------------
! read namelist
!-----------------------------------------------------------------------

  read (input_nml_file, nml=vert_turb_nml, iostat=io)
  ierr = check_nml_error(io,'vert_turb_nml')

!-----------------------------------------------------------------------
! write version number and namelist to logfile
!-----------------------------------------------------------------------

  call write_version_number (version, tag)
  if (mpp_pe() == mpp_root_pe())                                       &
                                   write (stdlog(), nml=vert_turb_nml)

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
  iz = Mgrid%iz

  allocate (uatm (ibd:ied, jbd:jed))
  allocate (vatm (ibd:ied, jbd:jed))
  allocate (v2   (ibd:ied, jbd:jed))
  allocate (v2atu(ibd:ied, jbd:jed))
  allocate (v2atv(ibd:ied, jbd:jed))
  allocate (cdm  (ibd:ied, jbd:jed))
  allocate (v2old(ibd:ied, jbd:jed))

  uatm = 0. ; vatm = 0.
  v2atu = 0.; v2atv = 0.
  v2 = 0.
  v2old = 0.

  allocate (logz1z0(jbd:jed))
  allocate (logz1zh(jbd:jed))
  allocate (     pp(jbd:jed))
  allocate (     vs(jbd:jed))
  allocate (     zz(jbd:jed))

  logz1z0 = log(zmref(jbd:jed)/z0)  ! lowest mass level
  logz1zh = log(zmref(jbd:jed)/zh)

  allocate (dzref(jbd:jed, kbd:ked))

  do k=kbd+1,ked-1
     dzref(:,k) = 0.5*(zref(jbd:jed,k-1) - zref(jbd:jed,k+1))
  enddo
  dzref(:,kbd) = zref(jbd:jed,kbd) - zref(jbd:jed,kbd+1)
  dzref(:,ked) = zref(jbd:jed,ked-1) - zref(jbd:jed,ked)

  do k=kbd,ked-1
     dzref(:,k) = dzref(:,k)/dzref(:,ked)
  enddo

  allocate (zetam(iz))
  allocate (zetaw(iz))
  allocate (dzetam(iz))
  allocate (dzetaw(iz))

  zetam = Vmetric%zetam
  zetaw = Vmetric%zetaw
  dzetam = Vmetric%dzetam
  dzetaw = Vmetric%dzetaw

  ptop = Vmetric%ptop

!-----------------------------------------------------------------------
! register diagnostics
!-----------------------------------------------------------------------

  axid_mm = (/axid(xm_),axid(ym_)/)

  name = "cdm"
  longname = "exchange coefficient for momentum"
  units = 'none'
  id_cdm = register_diag_field (                                       &
                       model, name, axid_mm, Time,                     &
                       longname, units, silly, range=(/-0.1, 0.1/) )
  name = "cdh"
  longname = "exchange coefficient for enthalpy"
  units = 'none'
  id_cdh = register_diag_field (                                       &
                       model, name, axid_mm, Time,                     &
                       longname, units, silly, range=(/-0.1, 0.1/) )

  return
end subroutine vert_turb_init

!#######################################################################

end module zetac_vert_turb_mod
