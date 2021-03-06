module zetac_user_ref_atmos_mod

!-----------------------------------------------------------------------
! module information
!-----------------------------------------------------------------------

use fms_mod,          only : check_nml_error, stdlog, stdout
use mpp_mod,          only : stdout, mpp_pe, input_nml_file

use constants_mod,    only : cp_air, rdgas, rvgas, hlv, grav

use zetac_convert_var_mod,       only : get_rh, get_pp
use zetac_coriolis_param_mod,    only : coriolis_param
use zetac_horiz_grid_type_mod,   only : horiz_grid_type
use zetac_horiz_metric_type_mod, only : horiz_metric_type
use zetac_moisture_mod,          only : get_esat, moisture_init
use zetac_vert_metric_type_mod,  only : vert_metric_type

!-----------------------------------------------------------------------
! public interfaces
!-----------------------------------------------------------------------

implicit none
public user_ref_atmos

character(len=*), parameter :: module='zetac_user_ref_atmos_mod'

!-----------------------------------------------------------------------
! version control information
!-----------------------------------------------------------------------

character(len=128)  :: version =  '$Id: zetac_user_ref_atmos.f90,v 1.1.2.9.2.15 2005/07/25 21:12:59 stg Exp $'
character(len=128)  :: tag     =  '$Name:  $'

contains

!#######################################################################

subroutine user_ref_atmos (                                            &
                                              Mgrid, Hmetric, Vmetric, &
                                    psref, tsref, uuref, taref, rhref )

!-----------------------------------------------------------------------
! program to  2D wind, temperature, pressure and moisture fields
!-----------------------------------------------------------------------

type (horiz_grid_type),                              intent (in)   ::  &
                                                                Mgrid

type (horiz_metric_type),                            intent(in)    ::  &
                                                              Hmetric

type (vert_metric_type),                             intent(in)    ::  &
                                                              Vmetric

real, dimension(Mgrid%jbg:Mgrid%jeg),                intent (out)  ::  &
                                                         psref, tsref

real, dimension(Mgrid%jbg:Mgrid%jeg, Mgrid%iz),      intent (out)  ::  &
                                                  uuref, taref, rhref

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

integer :: j, iy, jbg, jeg
integer :: k, iz
integer :: iymid
integer :: ierr, io

real, dimension(Mgrid%jbg:Mgrid%jeg, Mgrid%iz) ::                      &
                                                      corforce, ppref
real, dimension(Mgrid%jbg:Mgrid%jeg)           ::                      &
                                             ulat, vlat, ffv, cormetv

real :: dy
real :: expo, ptop, fac

!namelist:

real ::                     &
  lapse_rate = 0.,          &
  tsurf = 300.,             &
  diseq = 1.,               &
  rh0 = 0.7

namelist / reference_nml /  &
  lapse_rate,               &
  tsurf, diseq, rh0

real :: e0, p0, q0, t0

  jbg  = Mgrid%jbg
  jeg  = Mgrid%jeg
  iy   = Mgrid%iy
  iz   = Mgrid%iz
  
  ulat = Hmetric%ulat - Hmetric%rlatmin
  vlat = Hmetric%vlat - Hmetric%rlatmin
  dy   = Hmetric%dy

  iymid = (jbg + jeg)/2

!-----------------------------------------------------------------------
! read namelist
!-----------------------------------------------------------------------

  read (input_nml_file, nml=reference_nml, iostat=io)
  ierr = check_nml_error(io,'reference_nml')

!-----------------------------------------------------------------------
! get coriolis and spherical metric parameters
!-----------------------------------------------------------------------

  call coriolis_param ( Mgrid, Hmetric, vlat, ffv, cormetv )

!-----------------------------------------------------------------------
! get reference surface pressure
!-----------------------------------------------------------------------

  psref = Vmetric%pbot
  tsref = tsurf
  uuref = 0.
  rhref = rh0

!-----------------------------------------------------------------------
! get temperature profile at mass levels
!-----------------------------------------------------------------------

  ptop = Vmetric%ptop

  call get_pp ( psref, Vmetric%zetam, ptop, ppref )

  expo = Rdgas/Cp_air             ! dry adiabatic
  do j=jbg,jeg
     taref(j,iz) = (tsref(j) - diseq)*(ppref(j,iz)/psref(j))**expo
  enddo
  expo = Rdgas*lapse_rate/Grav    ! namelist lapse rate
  do k=iz-1,1,-1
     do j=jbg,jeg
        taref(j,k) = taref(j,k+1)*(ppref(j,k)/ppref(j,k+1))**expo
     enddo
  enddo
!!$if (mpp_pe() == 136) print *,psref(3), tsurf
if (mpp_pe() == 0) print *,"ppref",ppref(jbg,2:iz)
if (mpp_pe() == 0) print *,"taref",taref(jbg,2:iz)
!!$if (mpp_pe() == 136) print *,300.*(ppref(3,4:2:-1)/101300.)**(rdgas/cp_air)
!!$if (mpp_pe() == 136) print *,ppref(3,4:2:-1)

  taref(:,1) = taref(:,2)

  call moisture_init

  p0 = Vmetric%pbot

  t0 = tsurf - diseq
  call get_esat ( t0, e0 )
  q0 = (Rdgas/Rvgas)*e0/(p0 - e0)
  expo = q0

  t0 = 300. - diseq
  call get_esat ( t0, e0 )
  q0 = (Rdgas/Rvgas)*e0/(p0 - e0)

  diseq = rh0 * Hlv*(expo - q0)/Cp_air

  do k=1,iz
!     if (Vmetric%zetam(k) > .7) exit
     fac = 1. - dim(Vmetric%zetam(k), 0.2)/(Vmetric%zetam(iz) - 0.2)
     taref(:,k) = taref(:,k) + diseq*fac
  enddo
  if (mpp_pe() == 0) print *, "del T =", diseq

  return
end subroutine user_ref_atmos

!#######################################################################

end module zetac_user_ref_atmos_mod
