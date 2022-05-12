module zetac_convert_var_mod

!-----------------------------------------------------------------------
! module information
!-----------------------------------------------------------------------

use mpp_mod,                     only : mpp_pe, mpp_root_pe
use fms_mod,                     only : write_version_number
use constants_mod,               only : rdgas, rvgas, cp_air, kappa,   &
                                        hlv, tfreeze, grav
use zetac_moisture_mod,          only : get_qsat
use zetac_phys_con_mod,          only : rrat, eps, tnuke, tref, pref

!-----------------------------------------------------------------------
! public interfaces
!-----------------------------------------------------------------------

implicit none
private
public get_rh, get_qv, get_pp, get_tv, get_ro, get_sd, get_sm, get_phi

interface get_rh
  module procedure get_rh_3d
  module procedure get_rh_2d
  module procedure get_rh_0d
end interface

interface get_qv
  module procedure get_qv_3d
  module procedure get_qv_2d
end interface

interface get_pp
  module procedure get_pp_3d
  module procedure get_pp_2d
end interface

interface get_tv
  module procedure get_tv_3d
  module procedure get_tv_2d
end interface

interface get_ro
  module procedure get_ro_3d
  module procedure get_ro_2d
end interface

interface get_sd
  module procedure get_sd_3d
  module procedure get_sd_2d
end interface

interface get_sm
  module procedure get_sm_3d
  module procedure get_sm_2d
  module procedure get_sm_0d
end interface

interface get_phi
  module procedure get_phi_3d
  module procedure get_phi_2d
end interface

character(len=*), parameter :: module='zetac_convert_var_mod'
logical :: do_init=.true.

!-----------------------------------------------------------------------
! version control information
!-----------------------------------------------------------------------

character(len=128)  :: version =  '$Id: zetac_convert_var.f90,v 1.1.2.2.2.5 2005/07/05 02:01:37 stg Exp $'
character(len=128)  :: tag     =  '$Name:  $'

contains

!#######################################################################

subroutine get_rh_3d ( pp, ta, qv, rh, rhmin )

real, dimension(:,:,:),                     intent (in)    ::          &
                                                           pp, ta, qv

real, dimension(:,:,:),                     intent (out)   ::          &
                                                                   rh

real, optional,                             intent (in)    ::   rhmin

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

real, dimension(size(rh,1),size(rh,2),size(rh,3))          ::  qs, wi

!-----------------------------------------------------------------------
! initialize module
!-----------------------------------------------------------------------

  if (do_init) then
     call convert_var_init
  endif

!-----------------------------------------------------------------------
!  diagnose relative humidity
!-----------------------------------------------------------------------

  wi = tanh(dim(tfreeze-tnuke, ta)/tnuke)

!  call get_qsat ( pp, ta, qs, wi )
  call get_qsat ( pp, ta, qs )

  rh = min(1.1, qv/qs)

  if (present(rhmin)) rh = max(rh,rhmin)

  return
end subroutine get_rh_3d

!#######################################################################

subroutine get_rh_2d ( pp, ta, qv, rh, rhmin )

real, dimension(:,:),                       intent (in)    ::          &
                                                           pp, ta, qv

real, dimension(:,:),                       intent (out)   ::          &
                                                                   rh

real, optional,                             intent (in)    ::   rhmin

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

real, dimension(size(rh,1),size(rh,2))                     ::  qs, wi

!-----------------------------------------------------------------------
! initialize module
!-----------------------------------------------------------------------

  if (do_init) then
     call convert_var_init
  endif

!-----------------------------------------------------------------------
!  diagnose relative humidity
!-----------------------------------------------------------------------

  wi = tanh(dim(tfreeze-tnuke, ta)/tnuke)

!  call get_qsat ( pp, ta, qs, wi )
  call get_qsat ( pp, ta, qs )

  rh = qv/qs

  if (present(rhmin)) rh = max(rh,rhmin)

  return
end subroutine get_rh_2d

!#######################################################################

subroutine get_rh_0d ( pp, ta, qv, rh, rhmin )

real, intent (in)  :: pp, ta, qv
real, intent (out) :: rh

real, optional, intent (in) :: rhmin

real :: qs

  if (do_init) then
     call convert_var_init
  endif

  call get_qsat ( pp, ta, qs )

  rh = qv/qs

  if (present(rhmin)) rh = max(rh,rhmin)

  return
end subroutine get_rh_0d

!#######################################################################

subroutine get_qv_3d ( pp, ta, rh, qv )

real, dimension(:,:,:),                     intent (in)    ::          &
                                                           pp, ta, rh

real, dimension(:,:,:),                     intent (out)   ::          &
                                                                   qv

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------
 
real, dimension(size(qv,1),size(qv,2),size(qv,3))          ::  qs, wi

!-----------------------------------------------------------------------
! initialize module
!-----------------------------------------------------------------------

  if (do_init) then
     call convert_var_init
  endif

!-----------------------------------------------------------------------
!  diagnose specific humidity
!-----------------------------------------------------------------------

  wi = tanh(dim(tfreeze-tnuke, ta)/tnuke)  ! mixed-phase assumption

  call get_qsat ( pp, ta, qs, wi )

  qv = rh*qs
 
  return
end subroutine get_qv_3d

!#######################################################################

subroutine get_qv_2d ( pp, ta, rh, qv )

real, dimension(:,:),                       intent (in)    ::          &
                                                           pp, ta, rh

real, dimension(:,:),                       intent (out)   ::          &
                                                                   qv

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------
 
real, dimension(size(qv,1),size(qv,2))                     ::  qs, wi

!-----------------------------------------------------------------------
! initialize module
!-----------------------------------------------------------------------

  if (do_init) then
     call convert_var_init
  endif

!-----------------------------------------------------------------------
!  diagnose specific humidity
!-----------------------------------------------------------------------

  wi = tanh(dim(tfreeze-tnuke, ta)/tnuke)  ! mixed-phase assumption

  call get_qsat ( pp, ta, qs, wi )

  qv = rh*qs
 
  return
end subroutine get_qv_2d

!#######################################################################

subroutine get_pp_3d ( ps, zeta, ptop, pp )

real, dimension(:,:),                       intent (in)    ::          &
                                                                   ps

real, dimension(:),                         intent (in)    ::          &
                                                                 zeta

real,                                       intent (in)    ::          &
                                                                 ptop

real, dimension(:,:,:),                     intent (out)    ::         &
                                                                   pp

integer :: k

!-----------------------------------------------------------------------
! initialize module
!-----------------------------------------------------------------------

  if (do_init) then
     call convert_var_init
  endif

!-----------------------------------------------------------------------
! diagnose pressure from sigma
!-----------------------------------------------------------------------

  do k=1,size(zeta)
     pp(:,:,k) = ptop + (ps(:,:) - ptop)*zeta(k)
  enddo

  return
end subroutine get_pp_3d

!#######################################################################

subroutine get_pp_2d ( ps, zeta, ptop, pp )

real, dimension(:),                         intent (in)    ::          &
                                                                   ps

real, dimension(:),                         intent (in)    ::          &
                                                                 zeta

real,                                       intent (in)    ::          &
                                                                 ptop

real, dimension(:,:),                       intent (out)    ::         &
                                                                   pp

integer :: k

!-----------------------------------------------------------------------
! initialize module
!-----------------------------------------------------------------------

  if (do_init) then
     call convert_var_init
  endif

!-----------------------------------------------------------------------
! diagnose pressure from sigma
!-----------------------------------------------------------------------

  do k=1,size(zeta)
     pp(:,k) = ptop + (ps(:) - ptop)*zeta(k)
  enddo

  return
end subroutine get_pp_2d

!#######################################################################

subroutine get_tv_3d ( ta, qv, tv, qc )

real, dimension(:,:,:),                      intent (in)    ::         &
                                                               ta, qv

real, dimension(:,:,:),                      intent (out)   ::         &
                                                                   tv

real, optional, dimension(:,:,:),            intent (in)    ::         &
                                                                   qc

!-----------------------------------------------------------------------
! initialize module
!-----------------------------------------------------------------------

  if (do_init) then
     call convert_var_init
  endif

!-----------------------------------------------------------------------
! diagnose virtual temperature
!-----------------------------------------------------------------------

!  tv = ta*(1.0 + qv*rrat)
  tv = ta*(1.0 + qv/eps)/(1. + qv)

  if ( present (qc) ) tv = tv/(1.0 + qc)

  return
end subroutine get_tv_3d

!#######################################################################

subroutine get_tv_2d ( ta, qv, tv, qc )

real, dimension(:,:),                        intent (in)    ::         &
                                                               ta, qv

real, dimension(:,:),                        intent (out)   ::         &
                                                                   tv

real, optional, dimension(:,:),              intent (in)    ::         &
                                                                   qc

!-----------------------------------------------------------------------
! initialize module
!-----------------------------------------------------------------------

  if (do_init) then
     call convert_var_init
  endif

!-----------------------------------------------------------------------
! diagnose virtual temperature
!-----------------------------------------------------------------------

!  tv = ta*(1.0 + qv*rrat)
  tv = ta*(1.0 + qv/eps)/(1. + qv)

  if ( present (qc) ) tv = tv/(1.0 + qc)

  return
end subroutine get_tv_2d

!#######################################################################

subroutine get_ro_3d ( pp, tv, ro )

real, dimension(:,:,:),                      intent (in)    ::         &
                                                               pp, tv

real, dimension(:,:,:),                      intent (out)   ::         &
                                                                   ro

!-----------------------------------------------------------------------
! initialize module
!-----------------------------------------------------------------------

  if (do_init) then
     call convert_var_init
  endif

!-----------------------------------------------------------------------
! density from EOS
!-----------------------------------------------------------------------

  ro = pp/(Rdgas*tv)

end subroutine get_ro_3d

!#######################################################################

subroutine get_ro_2d ( pp, tv, ro )

real, dimension(:,:),                        intent (in)    ::         &
                                                               pp, tv

real, dimension(:,:),                        intent (out)   ::         &
                                                                   ro

!-----------------------------------------------------------------------
! initialize module
!-----------------------------------------------------------------------

  if (do_init) then
     call convert_var_init
  endif

!-----------------------------------------------------------------------
! density from EOS
!-----------------------------------------------------------------------

  ro = pp/(Rdgas*tv)

  return
end subroutine get_ro_2d

!#######################################################################

subroutine get_sd_3d ( ta, pp, sd )

real, dimension(:,:,:),                      intent (in)    ::         &
                                                               ta, pp

real, dimension(:,:,:),                      intent (out)   ::         &
                                                                   sd

!-----------------------------------------------------------------------
! initialize module
!-----------------------------------------------------------------------

  if (do_init) then
     call convert_var_init
  endif

!-----------------------------------------------------------------------
! moist entropy neglecting heat capacity
!-----------------------------------------------------------------------

  sd = Cp_air*log(ta/Tref) - Rdgas*log(pp/Pref)

end subroutine get_sd_3d

!#######################################################################

subroutine get_sd_2d ( ta, pp, sd )

real, dimension(:,:),                        intent (in)    ::         &
                                                               ta, pp

real, dimension(:,:),                        intent (out)   ::         &
                                                                   sd

!-----------------------------------------------------------------------
! initialize module
!-----------------------------------------------------------------------

  if (do_init) then
     call convert_var_init
  endif

!-----------------------------------------------------------------------
! dry entropy
!-----------------------------------------------------------------------

  sd = Cp_air*log(ta/Tref) - Rdgas*log(pp/Pref)

end subroutine get_sd_2d

!#######################################################################

subroutine get_sm_3d ( ta, pp, qv, rh, sm )

real, dimension(:,:,:),                      intent (in)    ::         &
                                                       ta, pp, qv, rh

real, dimension(:,:,:),                      intent (out)   ::         &
                                                                   sm

!-----------------------------------------------------------------------
! initialize module
!-----------------------------------------------------------------------

  if (do_init) then
     call convert_var_init
  endif

!-----------------------------------------------------------------------
! moist entropy neglecting heat capacity
!-----------------------------------------------------------------------

  sm = Cp_air*log(ta/Tref) - Rdgas*log(pp/(Pref*(1.+qv/eps)))          &
              + qv*(Hlv/ta - Rvgas*log(rh))

end subroutine get_sm_3d

!#######################################################################

subroutine get_sm_2d ( ta, pp, qv, rh, sm )

real, dimension(:,:),                        intent (in)    ::         &
                                                       ta, pp, qv, rh

real, dimension(:,:),                        intent (out)   ::         &
                                                                   sm

!-----------------------------------------------------------------------
! initialize module
!-----------------------------------------------------------------------

  if (do_init) then
     call convert_var_init
  endif

!-----------------------------------------------------------------------
! moist entropy neglecting heat capacity
!-----------------------------------------------------------------------

  sm = Cp_air*log(ta/Tref) - Rdgas*log(pp/(Pref*(1.+qv/eps)))          &
              + qv*(Hlv/ta - Rvgas*log(rh))

end subroutine get_sm_2d

!#######################################################################

subroutine get_sm_0d ( ta, pp, qv, rh, sm )

real, intent (in)  :: ta, pp, qv, rh

real, intent (out) :: sm

!-----------------------------------------------------------------------
! initialize module
!-----------------------------------------------------------------------

  if (do_init) then
     call convert_var_init
  endif

!-----------------------------------------------------------------------
! moist entropy neglecting heat capacity
!-----------------------------------------------------------------------

  sm = Cp_air*log(ta/Tref) - Rdgas*log(pp/(Pref*(1.+qv/eps)))          &
              + qv*(Hlv/ta - Rvgas*log(rh))

end subroutine get_sm_0d

!#######################################################################

subroutine get_phi_3d ( zsfc, pres, tvrt, geop )

!-----------------------------------------------------------------------
! get geopotential
!-----------------------------------------------------------------------

real, dimension(:,:),                                                  &
                                                intent (in)    ::      &
                                                                 zsfc

real, dimension(:,:,:),                                                &
                                                intent (in)    ::      &
                                                           pres, tvrt

real, dimension(:,:,:),                                                &
                                                intent (out)   ::      &
                                                                 geop

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

integer :: k, kdim

!-----------------------------------------------------------------------
! initialize module
!-----------------------------------------------------------------------

  if (do_init) then
     call convert_var_init
  endif

  kdim = size(pres,3)

  geop(:,:,kdim) = Grav*zsfc(:,:)

  do k=kdim,2,-1
     geop(:,:,k-1) = geop(:,:,k)                                       &
                    + Rdgas*tvrt(:,:,k)*log(pres(:,:,k)/pres(:,:,k-1))
  enddo

  return
end subroutine get_phi_3d

!#######################################################################

subroutine get_phi_2d ( pres, tvrt, geop )

!-----------------------------------------------------------------------
! get geopotential
!-----------------------------------------------------------------------

real, dimension(:,:),                                                  &
                                                intent (in)    ::      &
                                                           pres, tvrt

real, dimension(:,:),                                                  &
                                                intent (out)   ::      &
                                                                 geop

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

integer :: k, kdim

!-----------------------------------------------------------------------
! initialize module
!-----------------------------------------------------------------------

  if (do_init) then
     call convert_var_init
  endif

  kdim = size(pres,2)

  geop(:,kdim) = 0.

  do k=kdim,2,-1
     geop(:,k-1) = geop(:,k) + Rdgas*tvrt(:,k)*log(pres(:,k)/pres(:,k-1))
  enddo

  return
end subroutine get_phi_2d

!!#######################################################################

subroutine convert_var_init

!-----------------------------------------------------------------------
! write version number to logfile
!-----------------------------------------------------------------------

  call write_version_number (version, tag)

  do_init = .false.
  
  return
end subroutine convert_var_init

!#######################################################################

end module zetac_convert_var_mod

