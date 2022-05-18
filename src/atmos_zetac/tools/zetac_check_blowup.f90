module zetac_check_blowup_mod

!-----------------------------------------------------------------------
! module information
!-----------------------------------------------------------------------

use mpp_mod, only :   mpp_pe, mpp_root_pe, mpp_npes, mpp_sum,          &
                      input_nml_file
use fms_mod, only :   file_exist, write_version_number,                &
                      open_namelist_file, close_file,                  &
                      check_nml_error, stdlog, stdout
      
use zetac_horiz_grid_type_mod,   only : horiz_grid_type

!-----------------------------------------------------------------------
! public interfaces
!-----------------------------------------------------------------------

implicit none
private
public check_blowup

character(len=*), parameter :: module='zetac_check_blowup_mod'
logical :: do_init=.true.

real :: tolerance=0.0, tmin=180.0, tmax=340.0
real, allocatable, dimension(:,:) :: var1

integer :: off(3)
integer :: ibh, ieh
integer :: jbh, jeh

!-----------------------------------------------------------------------
! version control information
!-----------------------------------------------------------------------

character(len=128)  :: version =  '$Id: zetac_check_blowup.f90,v 1.1.2.8.2.2 2005/07/05 01:58:26 stg Exp $'
character(len=128)  :: tag     =  '$Name:  $'

contains

!#######################################################################

subroutine check_blowup ( Mgrid, var, temp, mxval, early_exit )

type (horiz_grid_type),                                                &
                                            intent(in)      ::  Mgrid

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),                                  &
                                            intent(in)      ::         &
                                                            var, temp

real,                                       intent(inout)   ::         &
                                                                mxval

logical,                                    intent(inout)   ::         &
                                                           early_exit

!-----------------------------------------------------------------------
!  local allocations
!-----------------------------------------------------------------------

integer :: loc(1)
real :: vmin, vmax, tcld, thot

!-----------------------------------------------------------------------
! initialize routine
!-----------------------------------------------------------------------

  if ( do_init ) then
     call check_blowup_init ( Mgrid )
  endif

  if ( tolerance <= 0.0 ) return

  var1 = 0.0

  var1(1,mpp_pe()) = minval(var(ibh:ieh,jbh:jeh,:))
  var1(2,mpp_pe()) = maxval(var(ibh:ieh,jbh:jeh,:))
  var1(3,mpp_pe()) = minval(temp(ibh:ieh,jbh:jeh,:))
  var1(4,mpp_pe()) = maxval(temp(ibh:ieh,jbh:jeh,:))

  call mpp_sum ( var1, 4*mpp_npes() )

  vmin = minval(var1(1,:))
  vmax = maxval(var1(2,:))
  tcld = minval(var1(3,:))
  thot = maxval(var1(4,:))

  mxval = sqrt(max(vmin*vmin, vmax*vmax))

  early_exit = ( mxval > tolerance .or. tcld < tmin .or. thot > tmax )

  if ( .not. early_exit ) return

  if ( mpp_pe() == 0) print *,"blowup: maxval =", mxval, tcld, thot
  loc = minloc(var1(1,:))
  if ( mpp_pe() == loc(1) ) then
     print *,"blowup: minval at", minloc(var(ibh:ieh,jbh:jeh,:)) + off
  endif
  loc = maxloc(var1(2,:))
  if ( mpp_pe() == loc(1) ) then
     print *,"blowup: maxval at", maxloc(var(ibh:ieh,jbh:jeh,:)) + off
  endif
  loc = minloc(var1(3,:))
  if ( mpp_pe() == loc(1) ) then
     print *,"coldest at", minloc(temp(ibh:ieh,jbh:jeh,:)) + off
  endif
  loc = maxloc(var1(4,:))
  if ( mpp_pe() == loc(1) ) then
     print *,"hottest at", maxloc(temp(ibh:ieh,jbh:jeh,:)) + off
  endif

  return
end subroutine check_blowup

!#######################################################################

subroutine check_blowup_init ( Mgrid )

type (horiz_grid_type),                                                &
                                            intent(in)      ::  Mgrid

!-----------------------------------------------------------------------
!  local allocations
!-----------------------------------------------------------------------

integer :: io, unit, ierr

namelist /blowup_nml/ tolerance, tmin, tmax

!-----------------------------------------------------------------------
! read namelist
!-----------------------------------------------------------------------

  read (input_nml_file, nml=blowup_nml, iostat=io)
  ierr = check_nml_error(io,'blowup_nml')

!-----------------------------------------------------------------------
! write version number and namelist to logfile
!-----------------------------------------------------------------------

  call write_version_number (version, tag)
  if (mpp_pe() == mpp_root_pe()) write (stdlog(), nml=blowup_nml)
  
  allocate (var1(4,0:mpp_npes()-1))

  ibh  = Mgrid%ibh
  ieh  = Mgrid%ieh
  jbh  = Mgrid%jbh
  jeh  = Mgrid%jeh

  off = (/ ibh-1 , jbh-1 , Mgrid%kbd-1 /)

  do_init = .false.

  return
end subroutine check_blowup_init

!#######################################################################

end module zetac_check_blowup_mod
