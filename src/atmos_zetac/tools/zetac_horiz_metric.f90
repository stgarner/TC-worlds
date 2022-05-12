module zetac_horiz_metric_mod

!-----------------------------------------------------------------------
! module information
!-----------------------------------------------------------------------

use mpp_mod,                     only : mpp_pe, mpp_root_pe
use fms_mod,                     only : write_version_number,          &
                                        error_mesg, FATAL
use constants_mod,               only : radian, radius

use zetac_horiz_grid_mod,        only : get_horiz_grid
use zetac_horiz_grid_type_mod,   only : horiz_grid_type
use zetac_horiz_metric_type_mod, only : horiz_metric_type,             &
                                        horiz_metric_type_define

!-----------------------------------------------------------------------
! public interfaces
!-----------------------------------------------------------------------

implicit none
private
public horiz_metric, cylindrical, cartesian, spherical

integer, parameter :: cylindrical=1, cartesian=2, spherical=3

character(len=*), parameter :: module='zetac_horiz_metric_mod'
logical :: do_init=.true.

real, parameter :: tiny=1.0e-12

!-----------------------------------------------------------------------
! version control information
!-----------------------------------------------------------------------

character(len=128)  :: version =  '$Id: zetac_horiz_metric.f90,v 1.1.2.2.2.6 2005/01/29 22:34:03 stg Exp $'
character(len=128)  :: tag     =  '$Name:  $'

contains

!#######################################################################

subroutine horiz_metric (                                              &
                                                       Mgrid, Hmetric, &
                                   rlonmin, rlonmax, rlatmin, rlatmax, &
                                         do_cartesian, do_cylindrical )

type (horiz_grid_type),                       intent (in)  ::          &
                                                                Mgrid

type (horiz_metric_type),                     intent (out) ::          &
                                                              Hmetric

real,                                         intent (in)  ::          &
                                   rlonmin, rlonmax, rlatmin, rlatmax
          
logical,                                      intent (in)  ::          &
                                         do_cartesian, do_cylindrical

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

real, dimension(Mgrid%ibg:Mgrid%ieg) :: ulon, vlon
real, dimension(Mgrid%jbg:Mgrid%jeg) :: ulat, vlat
real, dimension(Mgrid%jbg:Mgrid%jeg) :: cosulat, cosvlat
real, dimension(Mgrid%jbg:Mgrid%jeg) :: dxu, dxv
real, dimension(Mgrid%ibg:Mgrid%ieg, Mgrid%jbg:Mgrid%jeg) :: area

real    :: dlon, dlat, dx, dy, phi, rlat, dx1, areaj
integer :: geometry
integer :: i, ibg, ieg
integer :: j, jbg, jeg, iy

  ibg = Mgrid%ibg
  ieg = Mgrid%ieg
  jbg = Mgrid%jbg
  jeg = Mgrid%jeg

!-----------------------------------------------------------------------
! initialize module
!-----------------------------------------------------------------------

  if (do_init) then
     call horiz_metric_init
  endif

!-----------------------------------------------------------------------
! get horizontal grid 
!-----------------------------------------------------------------------

  call get_horiz_grid (                                                &
                                                                Mgrid, &
                                   rlonmin, rlonmax, rlatmin, rlatmax, &
                                                           dlon, dlat, &
                                               ulon, vlon, ulat, vlat )

  dx = dlon/radian*radius
  dy = dlat/radian*radius
  iy = Mgrid%iy

!-----------------------------------------------------------------------
! functions of latitude
!-----------------------------------------------------------------------

  if ( do_cylindrical ) then

     geometry = CYLINDRICAL

     if ( rlatmin < 0.0 ) call error_handler                           &
            ('rlatmin must be non-negative in cylindrical coordinates')

     do j=jbg,jeg
        cosulat(j) = max(tiny,ulat(j))/radian
        cosvlat(j) = max(tiny,vlat(j))/radian
        dxu(j) = dx*cosulat(j)
        dxv(j) = dx*cosvlat(j)
     enddo

  else if ( do_cartesian ) then

     geometry = CARTESIAN

     dx1 = dx*cos(0.5*(rlatmin + rlatmax)/radian)
     do j=jbg,jeg
        cosulat(j) = 1.0
        cosvlat(j) = 1.0
        dxu(j) = dx1
        dxv(j) = dx1
     enddo

  else

     geometry = SPHERICAL

     do j=jbg,jeg
        phi = ulat(j)/radian
        cosulat(j) = cos(phi)
        phi = vlat(j)/radian
        cosvlat(j) = cos(phi) 
        dxu(j) = dx*cosulat(j)
        dxv(j) = dx*cosvlat(j)
     enddo

     dx1 = dxu(1)
     where ( ulat(jbg:0) < -90.0 ) dxu(jbg:0) = dx1
     where ( vlat(jbg:0) < -90.0 ) dxv(jbg:0) = tiny
     dx1 = dxu(iy)
     where ( ulat(iy+1:jeg) > 90.0 ) dxu(iy+1:jeg) = dx1
     where ( vlat(iy+1:jeg) > 90.0 ) dxv(iy+1:jeg) = tiny

  endif

  dx = sum(dxu(1:iy))/iy

  do j=jbg,jeg
     areaj = dxu(j)*dy
     do i=ibg,ieg
        area(i,j) = areaj
     enddo
  enddo

!-----------------------------------------------------------------------
! define 'Hmetric'
!-----------------------------------------------------------------------

  call horiz_metric_type_define (                                      &
                                                       Mgrid, Hmetric, &
                                                             geometry, &
                                   rlonmin, rlonmax, rlatmin, rlatmax, &
                                               ulon, vlon, ulat, vlat, &
                                                     cosulat, cosvlat, &
                                         dlon, dlat, dx, dy, dxu, dxv, &
                                                                 area )

  return
end subroutine horiz_metric

!#######################################################################

subroutine horiz_metric_init

!-----------------------------------------------------------------------
! write version number to logfile
!-----------------------------------------------------------------------

  call write_version_number (version, tag)

  do_init = .false.
  
  return
end subroutine horiz_metric_init

!#######################################################################
   
subroutine error_handler ( message )

character(len=*), intent(in) :: message
   
  call error_mesg (module, message, FATAL)

  return
end subroutine error_handler

!#######################################################################

end module zetac_horiz_metric_mod

