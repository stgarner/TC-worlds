module zetac_axes_mod

!-----------------------------------------------------------------------
! module information
!-----------------------------------------------------------------------

use mpp_mod,          only : mpp_pe, mpp_root_pe
use mpp_io_mod,       only : axistype, mpp_write_meta, mpp_get_id
use fms_mod,          only : write_version_number
use mpp_domains_mod,  only : domain1D, mpp_get_domain_components,      &
                             mpp_get_global_domain, domain2D
use diag_manager_mod, only : diag_axis_init
use time_manager_mod, only : get_calendar_type, THIRTY_DAY_MONTHS,     &
                             JULIAN, GREGORIAN, NOLEAP

use zetac_axis_names_mod,        only : naxes, xu_, xv_, yu_, yv_,     &
                                               xm_, ym_, zw_, zm_,     &
                                               yg_, t_

use zetac_horiz_grid_type_mod,   only : horiz_grid_type
use zetac_horiz_metric_type_mod, only : horiz_metric_type
use zetac_vert_metric_type_mod,  only : vert_metric_type

!-----------------------------------------------------------------------
! public interfaces
!-----------------------------------------------------------------------

implicit none
private
public axes_history, axes_restart

character(len=*), parameter :: module='zetac_axes_mod'
!type(axistype), public, dimension(naxes), save :: Axes
integer, public, dimension(naxes) :: axid
integer, public, parameter :: gridu=1, gridv=2, gridm=3, gridw=4

logical :: do_init=.true.

type(domain1D), save :: Domain_x
type(domain1D), save :: Domain_y

!-----------------------------------------------------------------------
! version control information
!-----------------------------------------------------------------------

character(len=128)  :: version =  '$Id: zetac_axes.f90,v 1.1.2.4.2.4 2004/12/09 01:48:58 stg Exp $'
character(len=128)  :: tag     =  '$Name:  $'

contains

!#######################################################################

subroutine axes_history ( Mgrid, Hmetric, Vmetric )

type (horiz_grid_type),                     intent (in)    ::          &
                                                                Mgrid

type (horiz_metric_type),                   intent (in)    ::          &
                                                              Hmetric

type (vert_metric_type),                    intent (in)    ::          &
                                                              Vmetric

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

integer :: i, j, k, ix, iy
integer :: ibu, ieu, jbu, jeu
integer :: ibv, iev, jbv, jev

real, dimension(Mgrid%ibg:Mgrid%ieg) :: ulon, vlon
real, dimension(Mgrid%jbg:Mgrid%jeg) :: ulat, vlat
real, dimension(Mgrid%iz) :: zetaw
real, dimension(Mgrid%iz) :: zetam

  ix = Mgrid%ix
  iy = Mgrid%iy
  
  ulon = Hmetric%ulon
  vlon = Hmetric%vlon
  ulat = Hmetric%ulat
  vlat = Hmetric%vlat
  
  zetaw = Vmetric%zetaw
  zetam = Vmetric%zetam

!-----------------------------------------------------------------------
! initialize module
!-----------------------------------------------------------------------

  if (do_init) then
     call axes_init ( Mgrid )
  endif

!-----------------------------------------------------------------------
! set axis information
!-----------------------------------------------------------------------

  call mpp_get_global_domain (                                         &
                                                       Mgrid%Domain_u, &
                                                   ibu, ieu, jbu, jeu )

  call mpp_get_global_domain (                                         &
                                                       Mgrid%Domain_v, &
                                                   ibv, iev, jbv, jev )

  axid(xu_) = diag_axis_init ('lon_u', ulon(ibu:ieu), 'degrees_E', 'X',&
        'u_longitude', set_name='atmos', Domain2=Mgrid%Domain_u )

  axid(yu_) = diag_axis_init ('lat_u', ulat(jbu:jeu), 'degrees_N', 'Y',&
        'u_latitude' , set_name='atmos', Domain2=Mgrid%Domain_u )

  axid(xv_) = diag_axis_init ('lon_v', vlon(ibv:iev), 'degrees_E', 'X',&
        'v_longitude', set_name='atmos', Domain2=Mgrid%Domain_v )
 
  axid(yv_) = diag_axis_init ('lat_v', vlat(jbv:jev), 'degrees_N', 'Y',&
        'v_latitude' , set_name='atmos', Domain2=Mgrid%Domain_v )

  axid(xm_) = diag_axis_init ('lon_m', vlon(1:ix), 'degrees_E', 'X',   &
        'm_longitude', set_name='atmos', Domain2=Mgrid%Domain )
 
  axid(ym_) = diag_axis_init ('lat_m', ulat(1:iy), 'degrees_N', 'Y',   &
        'm_latitude' , set_name='atmos', Domain2=Mgrid%Domain )

  axid(yg_) = diag_axis_init ('lat_g', ulat(1:iy), 'degrees_N', 'Y',   &
        'global_lat', set_name='atmos' )

  axid(zw_) = diag_axis_init ('z_w', zetaw, 'none', 'Z',               &
        'w_level', set_name='atmos', direction=1 )

  axid(zm_) = diag_axis_init ('z_m', zetam, 'none', 'Z',               &
        'mass_level', set_name='atmos', direction=1 )

  return
end subroutine axes_history

!#######################################################################

subroutine axes_restart (                                              &
                                              Mgrid, Hmetric, Vmetric, &
                                                           unit, Axes, &
                                                     year, month, day )

type (horiz_grid_type),                    intent (in)    ::           &
                                                                Mgrid

type (horiz_metric_type),                  intent (in)    ::           &
                                                              Hmetric

type (vert_metric_type),                   intent (in)    ::           &
                                                              Vmetric

integer,                                   intent (in)    ::           &
                                                                 unit

type (axistype), dimension(:),             intent (inout)   ::           &
                                                                 Axes

integer, optional,                         intent (in)    ::           &
                                                     year, month, day

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

real, dimension(Mgrid%ibg:Mgrid%ieg) :: ulon, vlon
real, dimension(Mgrid%jbg:Mgrid%jeg) :: ulat, vlat
real, dimension(Mgrid%iz) :: zetaw
real, dimension(Mgrid%iz) :: zetam
character(len=2)  :: mstring, dstring
character(len=34) :: string
integer :: n, calendar, id
character(len=64) :: calendar_type, name, units
integer :: length

  ulon = Hmetric%ulon
  vlon = Hmetric%vlon
  ulat = Hmetric%ulat
  vlat = Hmetric%vlat

  zetaw = Vmetric%zetaw
  zetam = Vmetric%zetam

!-----------------------------------------------------------------------
! initialize module
!-----------------------------------------------------------------------

  if (do_init) then
     call axes_init ( Mgrid )
  endif

!-----------------------------------------------------------------------
! write axis information
!-----------------------------------------------------------------------

  do n=1,naxes

  if (n == xu_) then
     call mpp_write_meta( unit, Axes(xu_ ), 'lon_u', 'degrees_E',      &
                            'u_longitude', domain=Domain_x, data=ulon )

  else if (n == yu_) then
     call mpp_write_meta( unit, Axes(yu_), 'lat_u', 'degrees_N',       &
                            'u_latitude',  domain=Domain_y, data=ulat )

  else if (n == xv_) then
     call mpp_write_meta( unit, Axes(xv_ ), 'lon_v', 'degrees_E',      &
                            'v_longitude', domain=Domain_x, data=vlon )

  else if (n == yv_) then
     call mpp_write_meta( unit, Axes(yv_), 'lat_v', 'degrees_N',       &
                            'v_latitude',  domain=Domain_y, data=vlat )

  else if (n == xm_) then
     call mpp_write_meta( unit, Axes(xm_), 'lon_m', 'degrees_E',       &
                            'm_longitude', domain=Domain_x, data=vlon )

  else if (n == ym_) then
     call mpp_write_meta( unit, Axes(ym_), 'lat_m', 'degrees_N',       &
                            'm_latitude',  domain=Domain_y, data=ulat )

  else if (n == zw_) then
     call mpp_write_meta( unit, Axes(zw_), 'z_w', 'none',              &
                            'w_level',  data=zetaw )

  else if (n == zm_) then
     call mpp_write_meta( unit, Axes(zm_ ), 'z_m', 'none',             &
                            'mass_level', data=zetam )
if (mpp_pe() == 0) print *,zetam
  else if (n == yg_) then
     call mpp_write_meta( unit, Axes(yg_), 'lat', 'degrees_N',         &
                            'mass_latitude', data=ulat )

  else if (n == t_ .and. present(year) ) then
     if ( month <  10 ) write ( mstring, '("0",i1)' ) month
     if ( month >= 10 ) write ( mstring, '(i2)' ) month
     if ( day <  10 ) write ( dstring, '("0",i1)' ) day
     if ( day >= 10 ) write ( dstring, '(i2)' ) day
     write ( string, '("minutes since ",i4,"-",a2,"-",a2,              &
&                      " 00:00:00")' ) year, mstring, dstring
     call mpp_write_meta( unit, Axes(t_), 'time', string, 'time' )

     id = mpp_get_id ( Axes(t_) )
     calendar = get_calendar_type ( )

     select case (calendar)

        case (THIRTY_DAY_MONTHS)
           call mpp_write_meta ( unit, id,                             &
                             'calendar_type', cval='THIRTY_DAY_MONTHS')
           call mpp_write_meta ( unit, id,                             &
                             'calendar', cval='360')
        case (JULIAN)
           call mpp_write_meta ( unit, id,                             &
                             'calendar_type', cval='JULIAN')
           call mpp_write_meta ( unit, id,                             &
                             'calendar', cval='JULIAN')
        case (GREGORIAN)
           call mpp_write_meta ( unit, id,                             &
                             'calendar_type', cval='GREGORIAN')
           call mpp_write_meta ( unit, id,                             &
                             'calendar', cval='GREGORIAN')
        case (NOLEAP)
           call mpp_write_meta ( unit, id,                             &
                             'calendar_type', cval='NOLEAP')
           call mpp_write_meta ( unit, id,                             &
                             'calendar', cval='NOLEAP')

     end select

  endif

  enddo

  return
end subroutine axes_restart

!#######################################################################

subroutine axes_init ( Mgrid )

type (horiz_grid_type),                     intent (in)    ::          &
                                                                Mgrid

!-----------------------------------------------------------------------
! write version number to logfile
!-----------------------------------------------------------------------

  call write_version_number (version, tag)

  call mpp_get_domain_components ( Mgrid%Domain_plusxy,                &
                                                   Domain_x, Domain_y )

  do_init = .false.

  return
end subroutine axes_init

!#######################################################################

end module zetac_axes_mod
