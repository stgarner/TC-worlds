module zetac_user_perturb_mod

!-----------------------------------------------------------------------
! module information
!-----------------------------------------------------------------------

use mpp_mod,       only : mpp_pe, mpp_broadcast
use fms_mod,       only : file_exist, open_namelist_file, close_file

use zetac_horiz_grid_type_mod,    only : horiz_grid_type
use zetac_horiz_metric_type_mod,  only : horiz_metric_type
use zetac_rand_number_mod,        only : rand_number
use zetac_vert_metric_type_mod,   only : vert_metric_type

!-----------------------------------------------------------------------
! public interfaces
!-----------------------------------------------------------------------

implicit none

private
public user_perturb

character(len=*), parameter :: module='zetac_user_perturb_mod'
real :: seed=1.
!-----------------------------------------------------------------------
! version control information
!-----------------------------------------------------------------------

character(len=128)  :: version =  '$Id: zetac_user_perturb.f90,v 1.1.2.9.2.15 2005/07/25 21:12:59 stg Exp $'
character(len=128)  :: tag     =  '$Name:  $'

contains

!#######################################################################

subroutine user_perturb (                                              &
                                              Mgrid, Hmetric, Vmetric, &
                                                  psref, taref, rhref, &
                                    ppert, upert, vpert, tpert, hpert )

type (horiz_grid_type),   intent(in) ::                         Mgrid

type (horiz_metric_type), intent(in) ::                       Hmetric

type (vert_metric_type),  intent(in) ::                       Vmetric

real, dimension(Mgrid%jbg:Mgrid%jeg),                                  &
                                                    intent (in)   ::   &
                                                                psref

real, dimension(Mgrid%jbg:Mgrid%jeg,Mgrid%iz),                         &
                                                    intent (in)   ::   &
                                                         taref, rhref

real, dimension(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed),              &
                                                    intent (out)  ::   &
                                                                ppert

real, dimension(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed,               &
                Mgrid%kbd:Mgrid%ked),               intent (out)  ::   &
                                           upert, vpert, tpert, hpert

!-----------------------------------------------------------------------
!  local allocations
!-----------------------------------------------------------------------

integer :: iter, unit, io

integer :: i, ibh, ieh, ix
integer :: j, jbh, jeh, iy
integer :: k, kbd, ked
	   
real, dimension(Mgrid%iz) :: zetaw
real, dimension(Mgrid%ibd:Mgrid%ied) :: ulon
real, dimension(Mgrid%ibd:Mgrid%ied) :: vlon
real, dimension(Mgrid%jbd:Mgrid%jed) :: ulat
real, dimension(Mgrid%jbd:Mgrid%jed) :: vlat

real, dimension(Mgrid%iz) :: zeta

real ::                     &
  ample  = 0.0,             &
  xscale = 1.0e6,           &
  yscale = 1.0e6,           &
  zscale = 1.0

namelist / perturb_nml /    &
  ample, xscale, yscale, zscale

real :: xhome, yhome, zhome, xscal, yscal, zscal, facxy, facz

!-----------------------------------------------------------------------
!  initial non-zonal perturbation
!-----------------------------------------------------------------------

  ix  = Mgrid%ix
  iy  = Mgrid%iy

  ibh = Mgrid%ibh
  ieh = Mgrid%ieh
  jbh = Mgrid%jbh
  jeh = Mgrid%jeh

  kbd = Mgrid%kbd
  ked = Mgrid%ked

  ulat = Hmetric%ulat(Mgrid%jbd:Mgrid%jed)
  vlat = Hmetric%vlat(Mgrid%jbd:Mgrid%jed)
  ulon = Hmetric%ulon(Mgrid%ibd:Mgrid%ied)
  vlon = Hmetric%vlon(Mgrid%ibd:Mgrid%ied)

  zetaw = Vmetric%zetaw

!-----------------------------------------------------------------------
! read namelists
!-----------------------------------------------------------------------

  if( file_exist( 'input.nml' ) ) then
     unit = open_namelist_file ()
     io = 1
  do while ( io .ne. 0 )
     read( unit, nml = perturb_nml, iostat = io, end = 10 )
  end do
  10 call close_file (unit)
  endif

!-----------------------------------------------------------------------
! define 3D perturbation
!-----------------------------------------------------------------------

  ppert = 0.0
  upert = 0.0
  vpert = 0.0
  tpert = 0.0
  hpert = 0.0

  if (ample == 0.0) return

!-----------------------------------------------------------------------
!  pressure perturbation at mass levels between u and v points
!-----------------------------------------------------------------------

  xscal = xscale
  yscal = yscale
  zscal = zscale
  xhome = 0.5*(Hmetric%ulon(0) + Hmetric%ulon(ix))
  yhome = 0.5*(Hmetric%vlat(0) + Hmetric%vlat(iy))
  zhome = 0.5

  do k=kbd,ked
     facz = exp(-((zetaw(k)-zhome)/zscal)**2)
     do j=jbh,jeh
        do i=ibh,ieh
           facxy = exp(-((ulat(j) - yhome)/yscal)**2 - ((vlon(i) - xhome)/xscal)**2)
!           call rand_number(facxy)
           tpert(i,j,k) = ample*facxy*facz
        enddo
     enddo
  enddo

!-----------------------------------------------------------------------
! extrapolate perturbation at open lateral boundaries
!-----------------------------------------------------------------------

! call extrap_pert (                                                   &
!                                                               Mgrid, &
!                                   upert, vpert, epert, tpert, hpert )

  return
end subroutine user_perturb

!#######################################################################

end module zetac_user_perturb_mod
