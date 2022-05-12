module zetac_time_pointers_mod

!-----------------------------------------------------------------------
! module information
!-----------------------------------------------------------------------

use fms_mod, only : write_version_number

!-----------------------------------------------------------------------
! public interfaces
!-----------------------------------------------------------------------

implicit none
private
public update_time_pointers, get_time_pointers
public ntime

integer, parameter :: ntime=3

character(len=*), parameter :: module='zetac_time_pointers_mod'
logical :: do_init=.true.

integer :: past
integer :: present
integer :: future
integer :: sumppf

!-----------------------------------------------------------------------
! version control information
!-----------------------------------------------------------------------

character(len=128)  :: version =  '$Id: zetac_time_pointers.f90,v 1.1.2.1 2004/05/21 15:43:29 ck Exp $'
character(len=128)  :: tag     =  '$Name: latest $'

contains

!#######################################################################

subroutine update_time_pointers (   past_out, present_out, future_out, &
                                                           do_forward )

integer,                                            intent (out)  ::   &
                                    past_out, present_out, future_out

logical,                                            intent (in)   ::   &
                                                           do_forward
 
  past    = present
  present = future
  future  = sumppf - (past + present)

  if ( do_forward ) past = present

  past_out     = past
  present_out  = present
  future_out   = future

  return
end subroutine update_time_pointers

!#######################################################################
 
subroutine get_time_pointers ( past_get, present_get, future_get )

integer,                                           intent (out)   ::   &
                                     past_get, present_get, future_get

!-----------------------------------------------------------------------
! initialize time pointers if not done before
!-----------------------------------------------------------------------

  if ( do_init ) then
     call time_pointers_init
  endif
 
  past_get    = past
  present_get = present
  future_get  = future

  return
end subroutine get_time_pointers

!#######################################################################

subroutine time_pointers_init

!-----------------------------------------------------------------------
! write version number to logfile
!-----------------------------------------------------------------------

  call write_version_number (version, tag)

  past    = 1
  present = 2
  future  = 3

  sumppf = past + present + future
  
  do_init = .false.

  return
end subroutine time_pointers_init

!#######################################################################

end module zetac_time_pointers_mod

