module zetac_horiz_metric_type_mod

!-----------------------------------------------------------------------
! module information
!-----------------------------------------------------------------------

use fms_mod,                   only : write_version_number, mpp_pe

use zetac_horiz_grid_type_mod, only : horiz_grid_type

!-----------------------------------------------------------------------
! public interfaces
!-----------------------------------------------------------------------

implicit none
private
public horiz_metric_type, horiz_metric_type_define

type horiz_metric_type

  integer :: geometry

  real :: rlonmin, rlonmax, rlatmin, rlatmax
  real :: dx, dy, dlon, dlat

  real, pointer, dimension(:) :: ulon    =>NULL()
  real, pointer, dimension(:) :: vlon    =>NULL()
  real, pointer, dimension(:) :: ulat    =>NULL()
  real, pointer, dimension(:) :: vlat    =>NULL()
  real, pointer, dimension(:) :: cosulat =>NULL()
  real, pointer, dimension(:) :: cosvlat =>NULL()
  real, pointer, dimension(:) :: dxu     =>NULL()
  real, pointer, dimension(:) :: dxv     =>NULL()
  real, pointer, dimension(:,:) :: area  =>NULL()  !stg

end type horiz_metric_type

character(len=*), parameter :: module='zetac_horiz_metric_type_mod'
logical :: do_init=.true.

!-----------------------------------------------------------------------
! version control information
!-----------------------------------------------------------------------

character(len=128)  :: version =  '$Id: zetac_horiz_metric_type.f90,v 1.1.2.2.2.1 2004/06/15 04:15:18 stg Exp $'
character(len=128)  :: tag     =  '$Name: latest $'

contains

!#######################################################################

subroutine horiz_metric_type_define (                                  &
                                                       Mgrid, Hmetric, &
                                                             geometry, &
                                   rlonmin, rlonmax, rlatmin, rlatmax, &
                                               ulon, vlon, ulat, vlat, &
                                                     cosulat, cosvlat, &
                                         dlon, dlat, dx, dy, dxu, dxv, &
                                                                 area )

type (horiz_grid_type),               intent (in) ::            Mgrid

type (horiz_metric_type),             intent(out) ::          Hmetric

integer,                              intent (in) ::         geometry

real,                                 intent (in) :: rlonmin, rlonmax
real,                                 intent (in) :: rlatmin, rlatmax
real,                                 intent (in) ::       dlon, dlat
real,                                 intent (in) ::           dx, dy

real, dimension(Mgrid%ibg:Mgrid%ieg), intent (in) ::       ulon, vlon
real, dimension(Mgrid%jbg:Mgrid%jeg), intent (in) ::       ulat, vlat
real, dimension(Mgrid%jbg:Mgrid%jeg), intent (in) :: cosulat, cosvlat
real, dimension(Mgrid%jbg:Mgrid%jeg), intent (in) ::              dxu
real, dimension(Mgrid%jbg:Mgrid%jeg), intent (in) ::              dxv
real, dimension(Mgrid%ibg:Mgrid%ieg, Mgrid%jbg:Mgrid%jeg),             &
                                      intent (in) ::             area

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

  allocate (Hmetric%ulon   (Mgrid%ibg:Mgrid%ieg))
  allocate (Hmetric%vlon   (Mgrid%ibg:Mgrid%ieg))
  allocate (Hmetric%ulat   (Mgrid%jbg:Mgrid%jeg))
  allocate (Hmetric%vlat   (Mgrid%jbg:Mgrid%jeg))
  allocate (Hmetric%cosulat(Mgrid%jbg:Mgrid%jeg))
  allocate (Hmetric%cosvlat(Mgrid%jbg:Mgrid%jeg))
  allocate (Hmetric%dxu    (Mgrid%jbg:Mgrid%jeg))
  allocate (Hmetric%dxv    (Mgrid%jbg:Mgrid%jeg))
  allocate (Hmetric%area   (Mgrid%ibg:Mgrid%ieg, Mgrid%jbg:Mgrid%jeg))

!-----------------------------------------------------------------------
! initialize module
!-----------------------------------------------------------------------

  if (do_init) then
     call horiz_metric_type_init
  endif

  Hmetric%geometry = geometry

  Hmetric%rlonmin = rlonmin
  Hmetric%rlonmax = rlonmax
  Hmetric%rlatmin = rlatmin
  Hmetric%rlatmax = rlatmax

  Hmetric%dlon = dlon
  Hmetric%dlat = dlat
  Hmetric%dx   = dx
  Hmetric%dy   = dy

  Hmetric%ulon = ulon
  Hmetric%vlon = vlon

  Hmetric%ulat = ulat
  Hmetric%vlat = vlat

  Hmetric%cosulat = cosulat
  Hmetric%cosvlat = cosvlat

  Hmetric%dxu = dxu
  Hmetric%dxv = dxv

  Hmetric%area = area

  return
end subroutine horiz_metric_type_define

!#######################################################################

subroutine horiz_metric_type_init

!-----------------------------------------------------------------------
! write version number to logfile
!-----------------------------------------------------------------------

  call write_version_number (version, tag)

  do_init = .false.
  
  return
end subroutine horiz_metric_type_init

!#######################################################################

end module zetac_horiz_metric_type_mod
    

