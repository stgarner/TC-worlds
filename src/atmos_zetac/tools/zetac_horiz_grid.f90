module zetac_horiz_grid_mod

!-----------------------------------------------------------------------
! module information
!-----------------------------------------------------------------------

use mpp_mod,                   only : mpp_pe, mpp_root_pe
use fms_mod,                   only : write_version_number

use zetac_horiz_grid_type_mod, only : horiz_grid_type

!-----------------------------------------------------------------------
! public interfaces
!-----------------------------------------------------------------------

implicit none
private
public get_horiz_grid

character(len=*), parameter :: module='zetac_horiz_grid_mod'
logical :: do_init=.true.

!-----------------------------------------------------------------------
! version control information
!-----------------------------------------------------------------------

character(len=128)  :: version =  '$Id: zetac_horiz_grid.f90,v 1.1.2.3.2.1 2004/07/09 00:50:01 stg Exp $'
character(len=128)  :: tag     =  '$Name: latest $'

contains

!#######################################################################

subroutine get_horiz_grid (                                            &
                                                                Mgrid, &
                                   rlonmin, rlonmax, rlatmin, rlatmax, &
                                                           dlon, dlat, &
                                               ulon, vlon, ulat, vlat )
      
type (horiz_grid_type),                         intent (in)  ::        &
                                                                Mgrid

real,                                           intent (in)  ::        &
                                   rlonmin, rlonmax, rlatmin, rlatmax


real,                                           intent (out) ::        &
                                                           dlon, dlat
 
real, dimension(Mgrid%ibg:Mgrid%ieg),           intent (out) ::        &
                                                           ulon, vlon

real, dimension(Mgrid%jbg:Mgrid%jeg),           intent (out) ::        &
                                                           ulat, vlat

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

integer :: i, ibg, ieg, ix
integer :: j, jbg, jeg, iy

  ibg = Mgrid%ibg
  ieg = Mgrid%ieg
  jbg = Mgrid%jbg
  jeg = Mgrid%jeg

  ix = Mgrid%ix
  iy = Mgrid%iy

!-----------------------------------------------------------------------
! initialize module
!-----------------------------------------------------------------------

  if (do_init) then
     call horiz_grid_init
  endif

!-----------------------------------------------------------------------
! longitudes
!-----------------------------------------------------------------------

  dlon = (rlonmax - rlonmin)/ix
  do i=ibg,ieg
     ulon(i) = rlonmin + dlon* i
     vlon(i) = rlonmin + dlon*(i - 0.5)
  enddo

!-----------------------------------------------------------------------
! latitudes
!-----------------------------------------------------------------------

  dlat = (rlatmax - rlatmin)/iy
  do j=jbg,jeg
     ulat(j) = rlatmin + dlat*(j - 0.5)
     vlat(j) = rlatmin + dlat* j
  enddo

  return
end subroutine get_horiz_grid

!#######################################################################

subroutine horiz_grid_init

!-----------------------------------------------------------------------
! write version number to logfile
!-----------------------------------------------------------------------

  call write_version_number (version, tag)

  do_init = .false.
  
  return
end subroutine horiz_grid_init

!#######################################################################

end module zetac_horiz_grid_mod
