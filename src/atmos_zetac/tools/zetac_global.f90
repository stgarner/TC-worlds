module zetac_global_mod

!-----------------------------------------------------------------------
! module information
!-----------------------------------------------------------------------

use mpp_mod,                     only : mpp_pe, mpp_npes, mpp_sum
use fms_mod,                     only : write_version_number
use constants_mod,               only : grav
use zetac_horiz_grid_type_mod,   only : horiz_grid_type

!-----------------------------------------------------------------------
! public interfaces
!-----------------------------------------------------------------------

implicit none
private
public global

interface global
  module procedure global_3d
  module procedure global_2d
end interface

real, allocatable, dimension(:) :: var

logical :: do_init=.true.

character(len=*), parameter :: module='zetac_global_mod'

character(len=128)  :: version =  '$Id: zetac_global.f90,v 1.1.2.10.2.9 2005/08/07 00:01:45 stg Exp $'
character(len=128)  :: tag     =  '$Name:  $'

contains

!#######################################################################

subroutine global_3d ( Mgrid, field, fldsum, mass )

type (horiz_grid_type),                         intent (in)    ::      &
                                                                 Mgrid

real, dimension(:,:,:),                                                &
                                                intent (in)    ::      &
                                                           field, mass

real, intent(out) :: fldsum

!-----------------------------------------------------------------------
! initialize module
!-----------------------------------------------------------------------

  if (do_init) then
     call global_init
  endif

  var = 0.
  var(mpp_pe()) = sum(field*mass)/Grav
  call mpp_sum (var, size(var))
  fldsum = sum(var)/(Mgrid%ix*Mgrid%iy)
!  if (mpp_pe() == 0) print *,"check", fldsum

  end subroutine global_3d
 
!#######################################################################

subroutine global_2d ( Mgrid, field, fldsum, mass )

type (horiz_grid_type),                         intent (in)    ::      &
                                                                 Mgrid

real, dimension(:,:),                                                  &
                                                intent (in)    ::      &
                                                                 field

real, dimension(:,:), optional,                                        &
                                                intent (in)    ::      &
                                                                  mass

real, intent(out) :: fldsum

!-----------------------------------------------------------------------
! initialize module
!-----------------------------------------------------------------------

  if (do_init) then
     call global_init
  endif

  var = 0.
  if (present(mass)) then
     var(mpp_pe()) = sum(field*mass)/Grav
  else
     var(mpp_pe()) = sum(field)
  endif

  call mpp_sum (var, size(var))
  fldsum = sum(var)/(Mgrid%ix*Mgrid%iy)
!  if (mpp_pe() == 0) print *,"check", fldsum

  end subroutine global_2d
 
!!#######################################################################

subroutine global_init

!-----------------------------------------------------------------------
! write version number to logfile
!-----------------------------------------------------------------------

  call write_version_number (version, tag)

  allocate (var(0:mpp_npes()-1))

  do_init = .false.
  
  return
end subroutine global_init

!#######################################################################

end module zetac_global_mod

