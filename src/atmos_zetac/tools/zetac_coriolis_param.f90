module zetac_coriolis_param_mod

!-----------------------------------------------------------------------
! module information
!-----------------------------------------------------------------------

use mpp_mod,       only : mpp_pe, mpp_root_pe,                         &
                          input_nml_file
use fms_mod,       only : file_exist, open_namelist_file, close_file,  &
                          write_version_number, check_nml_error,       &
                          stdlog, stdout
use constants_mod, only : omega, radian, radius

use zetac_horiz_grid_type_mod,   only : horiz_grid_type
use zetac_horiz_metric_mod,      only : CYLINDRICAL, CARTESIAN,        &
                                        SPHERICAL
use zetac_horiz_metric_type_mod, only : horiz_metric_type

!-----------------------------------------------------------------------
! public interfaces
!-----------------------------------------------------------------------

implicit none
private
public coriolis_param

character(len=*), parameter :: module='zetac_coriolis_param_mod'
logical :: do_init=.true.

! namelist:
real :: fconst=0.0, beta=0.0

!-----------------------------------------------------------------------
! version control information
!-----------------------------------------------------------------------

character(len=128)  :: version =  '$Id: zetac_coriolis_param.f90,v 1.1.2.3.2.5 2005/06/16 18:12:02 stg Exp $'
character(len=128)  :: tag     =  '$Name:  $'

contains

!#######################################################################

subroutine coriolis_param (                                            &
                                                       Mgrid, Hmetric, &
                                                 rlat, ff, cormet, ee )

type (horiz_grid_type),                         intent (in)  ::        &  
                                                                Mgrid

type (horiz_metric_type),                       intent (in)  ::        &
                                                              Hmetric

real, dimension(Mgrid%jbg:Mgrid%jeg),           intent (in)  ::        &
                                                                 rlat

real, dimension(Mgrid%jbg:Mgrid%jeg),           intent (out) ::        &
                                                           ff, cormet
                                    
real, optional, dimension(Mgrid%jbg:Mgrid%jeg), intent (out) ::        &
                                                                   ee
                                    
!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

integer :: j, jbg, jeg
real    :: rlatmid, dlat
real    :: rlatmin, rlatmax, rlatj

real, dimension(Mgrid%jbg:Mgrid%jeg) :: phi

  jbg = Mgrid%jbg
  jeg = Mgrid%jeg

  rlatmin = Hmetric%rlatmin
  rlatmax = Hmetric%rlatmax
  dlat = Hmetric%dlat

  if (do_init) then
     call coriolis_param_init
  endif

  select case (Hmetric%geometry)

     case ( CYLINDRICAL )
     do j=jbg,jeg
        rlatj = max(0.5*dlat, rlat(j) - rlatmin)
        phi(j) = rlatj/radian
        ff(j) = fconst
        cormet(j) = -1.0/(phi(j)*radius)
     enddo
     if ( present(ee) ) ee = 0.0

     case ( CARTESIAN )
     rlatmid = 0.5*(rlatmin + rlatmax)
     do j=jbg,jeg
        ff(j) = fconst + beta*(rlat(j) - rlatmid)*radius/radian
        cormet(j) = 0.0
     enddo
     if ( present(ee) ) ee = 0.0

     case ( SPHERICAL )
     do j=jbg,jeg
        rlatj = max(-90.0 + 0.5*dlat, min(90.0 - 0.5*dlat, rlat(j)))
        phi(j) = rlatj/radian
        ff(j) = 2.0*omega*sin(phi(j))
        cormet(j) = tan(phi(j))/radius
     enddo
     if ( present(ee) ) then
        do j=jbg,jeg
           ee(j) = 2.0*omega*cos(phi(j))
        enddo
     endif

  end select
  
  return
end subroutine coriolis_param

!#######################################################################

subroutine coriolis_param_init

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

integer :: io, unit, ierr

namelist /coriolis_nml/ fconst, beta

!-----------------------------------------------------------------------
! read namelist
!-----------------------------------------------------------------------

  read (input_nml_file, nml=coriolis_nml, iostat=io)
  ierr = check_nml_error(io,'coriolis_nml')
 
!-----------------------------------------------------------------------
! write version number and namelist to logfile
!-----------------------------------------------------------------------

  call write_version_number (version, tag)
  if (mpp_pe() == mpp_root_pe())                                       &
                                     write (stdlog(), nml=coriolis_nml)
 
  do_init = .false.

  return
end subroutine coriolis_param_init

!#######################################################################

end module zetac_coriolis_param_mod
