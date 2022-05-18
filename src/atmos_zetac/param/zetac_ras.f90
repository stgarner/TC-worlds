module zetac_ras_mod

use zetac_axes_mod,         only : axid, gridm, gridw
use zetac_axis_names_mod,   only : xm_, ym_, zm_
use zetac_horiz_grid_type_mod, only : horiz_grid_type
use zetac_moisture_mod,     only : get_qsat, get_dlnesatdt
use zetac_ncdf_io_mod,      only : ncdf_fms_write
use zetac_phys_con_mod,     only : eps ! Rd/Rv
use fms_mod,                only : mpp_pe, mpp_root_pe, file_exist,    &
                                   error_mesg, check_nml_error, FATAL, &
                                   open_file, close_file,              &
                                   open_namelist_file, stdlog,         &
                                   write_version_number
use mpp_mod,                only : input_nml_file
use constants_mod,          only : cp_air, rdgas, rvgas, kappa, hlv,   &
                                   grav
use diag_manager_mod,       only : register_diag_field
use field_manager_mod,      only : model_atmos
use time_manager_mod,       only : time_type

implicit none
private

public ras, ras_init, ras_end

real,    parameter :: silly=-9999.9
real,    parameter :: lapse = grav/cp_air
integer, parameter :: max_nlam = 10

character (len=*),  parameter :: module='zetac_ras_mod'
character (len=16), parameter :: model='atmos_mod'

integer :: id_rain_conv, id_mc_conv
integer :: id_work_conv, id_tdt_conv, id_qdt_conv, id_ldt_conv
integer :: nlam

! namelist:
real :: evfac = 0.2
real :: work_crit(max_nlam)
real :: eff_ras(max_nlam)
real :: dt_ras(max_nlam)
real :: frac_ras(max_nlam)
logical :: use_ras = .true.

namelist /ras_nml/ evfac,                                              &
                   work_crit, dt_ras, eff_ras,                         &
                   use_ras

real, parameter :: locp = hlv/cp_air

real, parameter :: missing_value = -999.
real, parameter :: tiny = 1.e-10

character(len=16) :: mod_name = 'zetac_ras_mod'

character(len=128) :: version = '$Id: zetac_ras.f90,v 1.1.2.1.2.2.2.2 2004/09/15 16:03:45 ck Exp $'
character(len=128) :: tagname = '$Name:  $'

contains

!#######################################################################

subroutine ras (                                                       &
                                                          Mgrid, Time, &
                                           zfull, pfull, zhalf, phalf, &
                                                           ta, qv, ql, &
                                                                 delt )
 
type (horiz_grid_type),                      intent (in)    ::         &
                                                                Mgrid
       
type (time_type),                            intent (in)    ::         &
                                                                 Time

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),        intent(in)     ::         &
                                           zfull, pfull, zhalf, phalf
  
real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                 Mgrid%kbd:Mgrid%ked),       intent(inout)  ::         &
                                                           ta, qv, ql
 

real,                                        intent(in)     ::   delt

!-----------------------------------------------------------------------
! local allocations  
!-----------------------------------------------------------------------

integer :: l, k, iz

real, dimension(size(ta,1),size(ta,2),size(ta,3)) :: lam, eta, work, hc, qc
real, dimension(size(ta,1),size(ta,2),size(ta,3)) :: tdt, qdt, ldt, gam
real, dimension(size(ta,1),size(ta,2),size(ta,3)) :: ss, hh, hs, hint, qint
real, dimension(size(ta,1),size(ta,2),size(ta,3)) :: qsat, dqdt
real, dimension(size(ta,1),size(ta,2),size(ta,3)) :: dsdp, dhdp
real, dimension(size(ta,1),size(ta,2)) :: kern, mflux, mc
real, dimension(size(ta,1),size(ta,2)) :: qcon, rain

  iz = size(ta,3)

  call get_qsat ( pfull, ta, qsat )
  call get_dlnesatdt ( ta, dqdt )
  dqdt = qsat * dqdt * (1. + qsat/eps)
  gam(:,:,:) = 1.0 + locp*dqdt(:,:,:)

  call get_lambda ( zfull, ta, qv, qsat, ss, hh, hs, hint, qint, lam )
  call get_stab ( pfull, ss, hh, dsdp, dhdp )

!if (mpp_pe() == 137) print *, "lam",lam(3,3,iz-1:2:-1)

  mc = 0.
  work = 0.
  tdt = 0.
  qdt = 0.
  rain = 0.
ldt=0.

  do l=iz-1,2,-1

! mass-flux profile

     do k=iz,l,-1
        eta(:,:,k) = 1.0 + lam(:,:,l)*(zfull(:,:,k) - zfull(:,:,iz))
     enddo

! workfunction (units: m2/s2)

     call workfunction ( l, lam(:,:,l), zfull, pfull, gam, eta,       &
                                hs, hint, qint, hh, qv, hc, qc, work )

! kernel (units: m2/s2/Pa) and mass flux (units: Pa/s)

     call kernel ( l, lam(:,:,l), zfull, pfull,                       &
                          qc, dsdp, dhdp, eta, gam, qsat, qcon, kern )

     where ( kern(:,:) < 0.1 .or. lam(:,:,l) > 0.001                  &
                             .or. eta(:,:,l) < 0) work(:,:,l) = 0.

     mflux(:,:) = dim(work(:,:,l),work_crit(l))/(kern(:,:)*dt_ras(l))

! tendencies

     do k=iz,l,-1
        qdt(:,:,k) = qdt(:,:,k) + (dhdp(:,:,k) - dsdp(:,:,k))*mflux(:,:)/hlv
        tdt(:,:,k) = tdt(:,:,k) + dsdp(:,:,k)*mflux(:,:)/cp_air
     enddo

! handle detrained condensate

     qcon(:,:) = qcon(:,:)*mflux(:,:)
     ldt(:,:,l) = qcon(:,:)/(phalf(:,:,l) - phalf(:,:,l-1))*          &
                                                    (1.0 - eff_ras(l))

! accumulate convective mass flux and rain

     mc(:,:) = mc(:,:) + mflux(:,:)
     rain(:,:) = rain(:,:) + qcon(:,:)/Grav*eff_ras(l)

  enddo

  if (use_ras) then
     ta(:,:,:) = ta(:,:,:) + tdt(:,:,:)*delt
     qv(:,:,:) = qv(:,:,:) + qdt(:,:,:)*delt
     ql(:,:,:)=  ql(:,:,:) + ldt(:,:,:)*delt
  endif

!-----------------------------------------------------------------------
! write diagnostic fields
!-----------------------------------------------------------------------

  call ncdf_fms_write (Mgrid, Time, id_rain_conv, gridm, rain)
  call ncdf_fms_write (Mgrid, Time, id_mc_conv, gridm, mc)
  call ncdf_fms_write (Mgrid, Time, id_work_conv, gridm, work)
  call ncdf_fms_write (Mgrid, Time, id_tdt_conv, gridm, tdt)
  call ncdf_fms_write (Mgrid, Time, id_qdt_conv, gridm, qdt)
  call ncdf_fms_write (Mgrid, Time, id_ldt_conv, gridm, ldt)

  return
 end subroutine ras

!#######################################################################

 subroutine get_lambda ( zfull, ta, qv, qsat, ss, hh, hs, hint, qint, lam )

  real, intent(in),  dimension(:,:,:) :: zfull
  real, intent(in),  dimension(:,:,:) :: ta, qv, qsat
  real, intent(out), dimension(:,:,:) :: ss, hh, hs, hint, qint, lam

  integer :: k, iz

!-----------------------------------------------------------------------
! entrainment parameter for given detrainment level (units: 1/length)
!-----------------------------------------------------------------------

  iz = size(ta,3)

  do k=2,iz
     ss(:,:,k) = cp_air*ta(:,:,k) + grav*zfull(:,:,k)
     hh(:,:,k) = ss(:,:,k) + hlv*qv(:,:,k)
     hs(:,:,k) = ss(:,:,k) + hlv*qsat(:,:,k)
  enddo

  hint(:,:,iz) = 0.
  qint(:,:,iz) = 0.

  do k=iz-1,2,-1
     hint(:,:,k) = hint(:,:,k+1) + 0.5*(hh(:,:,k+1) + hh(:,:,k))*      &
                                    (zfull(:,:,k) - zfull(:,:,k+1))
     qint(:,:,k) = qint(:,:,k+1) + 0.5*(qv(:,:,k+1) + qv(:,:,k))*      &
                                    (zfull(:,:,k) - zfull(:,:,k+1))
  enddo

  lam(:,:,1) = 0.
  lam(:,:,iz) = 0.

  do k=2,iz-1
     lam(:,:,k) = dim(hh(:,:,iz), hs(:,:,k)) / max( tiny,              &  ! Eq. 10
              hs(:,:,k)*(zfull(:,:,k) - zfull(:,:,iz)) - hint(:,:,k) )
  enddo

 end subroutine get_lambda 

!#######################################################################

 subroutine get_stab ( pfull, ss, hh, dsdp, dhdp )

  real, intent(in),  dimension(:,:,:) :: pfull
  real, intent(in),  dimension(:,:,:) :: ss, hh
  real, intent(out), dimension(:,:,:) :: dsdp, dhdp

  integer :: k, iz

  iz = size(hh,3)

!-----------------------------------------------------------------------
! stability (units: m2/s2/Pa)
!-----------------------------------------------------------------------

  do k=3,iz-1
     dsdp(:,:,k) =                                                     &
          (ss(:,:,k-1) - ss(:,:,k+1))/(pfull(:,:,k+1) - pfull(:,:,k-1))
     dhdp(:,:,k) =                                                     &
          (hh(:,:,k-1) - hh(:,:,k+1))/(pfull(:,:,k+1) - pfull(:,:,k-1))
  enddo

  dsdp(:,:,2) = 1.0*(ss(:,:,2) - ss(:,:,3))/                           &
                    (pfull(:,:,3) - pfull(:,:,2))
  dhdp(:,:,2) = 1.0*(hh(:,:,2) - hh(:,:,3))/                           &
                    (pfull(:,:,3) - pfull(:,:,2))
  dsdp(:,:,iz) = 0.5*(ss(:,:,iz-1) - ss(:,:,iz))/                      &
                    (pfull(:,:,iz) - pfull(:,:,iz-1))
  dhdp(:,:,iz) = 0.5*(hh(:,:,iz-1) - hh(:,:,iz))/                      &
                    (pfull(:,:,iz) - pfull(:,:,iz-1))

end subroutine get_stab

!#######################################################################

subroutine workfunction ( level, lam, zfull, pfull, gam, eta,          &
                                 hs, hint, qint, hh, qv, hc, qc, work )

integer, intent(in) :: level
real, intent(in),  dimension(:,:)   :: lam
real, intent(in),  dimension(:,:,:) :: zfull, pfull, gam, eta
real, intent(in),  dimension(:,:,:) :: hs, hint, qint, hh, qv
real, intent(out), dimension(:,:,:) :: hc, qc, work

integer :: k, iz
real, dimension(size(hh,1),size(hh,2),size(hh,3)) :: buoy

  iz = size(hs,3)

!-----------------------------------------------------------------------
! compute cloud work function (units: m2/s2)
!-----------------------------------------------------------------------

  do k=iz,level,-1
     hc(:,:,k) = (hh(:,:,iz) + lam(:,:)*hint(:,:,k))/eta(:,:,k)            ! Eq. 9
     qc(:,:,k) = (qv(:,:,iz) + lam(:,:)*qint(:,:,k))/eta(:,:,k)
  enddo

  do k=iz,level,-1
     buoy(:,:,k) = eta(:,:,k)*(hc(:,:,k) - hs(:,:,k))/gam(:,:,k)
  enddo
  work(:,:,level) = 0.
  do k=iz-1,level,-1
     work(:,:,level) = work(:,:,level)                                 &  !  Eq. 13
                                      + (buoy(:,:,k+1) + buoy(:,:,k))* &
                                      (pfull(:,:,k+1) - pfull(:,:,k))/ &
                                      (pfull(:,:,k+1) + pfull(:,:,k))
  enddo

end subroutine workfunction

!#######################################################################

 subroutine kernel ( level, lam, zfull, pfull,                         &
                           qc, dsdp, dhdp, eta, gam, qsat, qcon, kern )

  integer, intent(in) :: level
  real, intent(in),    dimension(:,:) :: lam
  real, intent(in),  dimension(:,:,:) :: zfull, pfull, qc
  real, intent(in),  dimension(:,:,:) :: dsdp, dhdp, eta, gam
  real, intent(in),  dimension(:,:,:) :: qsat
  real, intent(out),   dimension(:,:) :: qcon, kern

  integer :: k, iz
  real, dimension(size(qc,1),size(qc,2),size(qc,3)) :: hhtend
  real, dimension(size(qc,1),size(qc,2),size(qc,3)) :: hstend
  real, dimension(size(qc,1),size(qc,2),size(qc,3)) :: htint
  real, dimension(size(qc,1),size(qc,2),size(qc,3)) :: buoy

  iz = size(qc,3)

!-----------------------------------------------------------------------
! compute mass-flux kernel (units: m2/s2/Pa)
!-----------------------------------------------------------------------

  do k=iz,level,-1
     hhtend(:,:,k) = eta(:,:,k)*dhdp(:,:,k)   ! tendency at full levels
     hstend(:,:,k) = eta(:,:,k)*dsdp(:,:,k)*gam(:,:,k)
  enddo

!-----------------------------------------------------------------------
! re-evaporate some cloud water at detrainment level
!-----------------------------------------------------------------------

  qcon(:,:) = dim(qc(:,:,level), qsat(:,:,level))*eta(:,:,level)

!-----------------------------------------------------------------------
! compute kernel
!-----------------------------------------------------------------------

  htint(:,:,iz) = 0.

  do k=iz-1,level,-1
     htint(:,:,k) = htint(:,:,k+1)                                     &
                              + 0.5*(hhtend(:,:,k) + hhtend(:,:,k+1))* &
                                      (zfull(:,:,k) - zfull(:,:,k+1))
  enddo

  do k=iz,level,-1
     buoy(:,:,k) = -(hhtend(:,:,iz) + lam(:,:)*htint(:,:,k)            &   
                   - hstend(:,:,k)*eta(:,:,k))/gam(:,:,k)
  enddo

  kern = 0.
  do k=iz-1,level,-1
     kern(:,:) = kern(:,:) + (buoy(:,:,k+1) + buoy(:,:,k))*            &
                            (pfull(:,:,k+1) - pfull(:,:,k))/           &
                            (pfull(:,:,k+1) + pfull(:,:,k))
  enddo

end subroutine kernel

!#######################################################################

subroutine ras_init ( Mgrid, Time )

type (horiz_grid_type),                     intent (in)    ::          &
                                                                Mgrid
       
type (time_type),                           intent (in)    ::          &
                                                                 Time

integer :: ierr, io
integer, dimension(2) :: axid_mm
integer, dimension(3) :: axid_mmm

character (len=16) :: name, units
character (len=64) :: longname

  dt_ras = 1800.
  work_crit = 0.
  frac_ras = 0.

!-----------------------------------------------------------------------
! namelist (read & write) 
!-----------------------------------------------------------------------

  read (input_nml_file, nml=ras_nml, iostat=io)
  ierr = check_nml_error(io,'ras_nml')

!-----------------------------------------------------------------------
! write version number and namelist to logfile
!-----------------------------------------------------------------------

  call write_version_number (version, tagname)
  if (mpp_pe() == mpp_root_pe()) write (stdlog(), nml=ras_nml)

  nlam = Mgrid%iz - 2

!-----------------------------------------------------------------------
! register diagnostics
!-----------------------------------------------------------------------

  axid_mmm = (/axid(xm_),axid(ym_),axid(zm_)/)
  axid_mm  = (/axid(xm_),axid(ym_)/)

  name = "rain_conv"
  longname = "rain rate from RAS"
  units = 'kg/m2/s'
  id_rain_conv = register_diag_field (                                 &
                      model, name, axid_mm, Time,                      &
                      longname, units, silly, range=(/-1.e-2, 1.e-2/) )

  name = "mc_conv"
  longname = "cloud mass flux from RAS"
  units = 'kg/m2/s'
  id_mc_conv = register_diag_field (                                   &
                      model, name, axid_mm, Time,                      &
                      longname, units, silly, range=(/-10., 10./) )

  name = "work_conv"
  longname = "work function from RAS"
  units = 'm2/s2'
  id_work_conv = register_diag_field (                                 &
                      model, name, axid_mmm, Time,                     &
                      longname, units, silly, range=(/-1.e6, 1.e6/) )

  name = "tdt_conv"
  longname = "temperature tendency from RAS"
  units = 'K/s'
  id_tdt_conv = register_diag_field (                                  &
                      model, name, axid_mmm, Time,                     &
                      longname, units, silly, range=(/-5.e-3, 5.e-3/) )

  name = "qdt_conv"
  longname = "specific humidity tendency from RAS"
  units = 'kg/kg/s'
  id_qdt_conv = register_diag_field (                                  &
                      model, name, axid_mmm, Time,                     &
                      longname, units, silly, range=(/-1.e-5, 1.e-5/) )

  name = "ldt_conv"
  longname = "condensate tendency from RAS"
  units = 'kg/kg/s'
  id_ldt_conv = register_diag_field (                                  &
                      model, name, axid_mmm, Time,                     &
                      longname, units, silly, range=(/-1.e-5, 1.e-5/) )

  return
end subroutine ras_init

!#######################################################################

subroutine ras_end ( Mgrid, time, ltime )

type (horiz_grid_type)                                       ::        &
                                                                Mgrid

real,                                         intent (in)    ::        &
                                                                 time

integer,                                      intent (in)    ::        &
                                                                ltime

  return
end subroutine ras_end

!#######################################################################

end module zetac_ras_mod

