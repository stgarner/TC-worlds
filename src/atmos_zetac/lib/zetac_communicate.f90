module zetac_communicate_mod

use mpp_mod, only : mpp_pe, mpp_sum, mpp_root_pe,                      &
                    mpp_send, mpp_recv, mpp_sync_self

!-----------------------------------------------------------------------
! public interface
!-----------------------------------------------------------------------

implicit none
private
public collect

interface collect
  module procedure collect_1d
  module procedure collect_2d
  module procedure collect_3d
  module procedure collect_4d
end interface

!-----------------------------------------------------------------------
! version control information
!-----------------------------------------------------------------------

character(len=128)  :: version =  '$Id: zetac_communicate.f90,v 1.1.2.4.2.2 2004/07/09 00:47:18 stg Exp $'
character(len=128)  :: tag     =  '$Name: latest $'

contains

!#######################################################################

subroutine collect_1d ( input, output, pelist )

real   , dimension(:), intent (inout) :: input
real   , dimension(:), intent (inout) :: output
integer, dimension(:), intent (in)    :: pelist

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

  integer :: n, pe, npes, endpe
  real, dimension(size(output,1)) :: temp

  pe    = mpp_pe()
  npes  = size(pelist)
  endpe = pelist(1)

  if ( pe /= endpe ) then
    call mpp_send ( input, size(input), endpe )
  else
    output = input
    do n=2,npes
      call mpp_recv ( input, size(input), pelist(n) )
      output = output + input
    enddo
  endif

  if ( pe == endpe ) then
    temp = output
    do n=2,npes
      call mpp_send ( temp, size(temp), pelist(n) )
    enddo
  else
    call mpp_recv ( output, size(output), endpe )
  endif

  call mpp_sync_self ( pelist )

  return
end subroutine collect_1d

!#######################################################################

subroutine collect_2d ( input, output, pelist )

real   , dimension(:,:), intent (inout) :: input
real   , dimension(:,:), intent (inout) :: output
integer, dimension(:)  , intent (in)    :: pelist

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

  integer :: n, pe, npes, endpe
  real, dimension(size(output,1),size(output,2)) :: temp

  pe    = mpp_pe()
  npes  = size(pelist)
  endpe = pelist(1)

  if ( pe /= endpe ) then
    call mpp_send ( input, size(input), endpe )
  else
    output = input
    do n=2,npes
      call mpp_recv ( input, size(input), pelist(n) )
      output = output + input
    enddo
  endif

  if ( pe == endpe ) then
    temp = output
    do n=2,npes
      call mpp_send ( temp, size(temp), pelist(n) )
    enddo
  else
    call mpp_recv ( output, size(output), endpe )
  endif

  call mpp_sync_self ( pelist )

  return
end subroutine collect_2d

!#######################################################################

subroutine collect_3d ( input, output, pelist )

real   , dimension(:,:,:), intent (inout) :: input
real   , dimension(:,:,:), intent (inout) :: output
integer, dimension(:)    , intent (in)    :: pelist

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

  integer :: n, pe, npes, endpe
  real, dimension(size(output,1),size(output,2),size(output,3)) :: temp

  pe    = mpp_pe()
  npes  = size(pelist)
  endpe = pelist(1)

  if ( pe /= endpe ) then
    call mpp_send ( input, size(input), endpe )
  else
    output = input
    do n=2,npes
      call mpp_recv ( input, size(input), pelist(n) )
      output = output + input
    enddo
  endif

  if ( pe == endpe ) then
    temp = output
    do n=2,npes
      call mpp_send ( temp, size(temp), pelist(n) )
    enddo
  else
    call mpp_recv ( output, size(output), endpe )
  endif

  call mpp_sync_self ( pelist )

  return
end subroutine collect_3d

!#######################################################################

subroutine collect_4d ( input, output, pelist )

real   , dimension(:,:,:,:), intent (inout) :: input
real   , dimension(:,:,:,:), intent (inout) :: output
integer, dimension(:)      , intent (in)    :: pelist

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

  integer :: n, pe, npes, endpe
  real, dimension(size(output,1),size(output,2),                       &
                  size(output,3),size(output,4)) :: temp

  pe    = mpp_pe()
  npes  = size(pelist)
  endpe = pelist(1)

  if ( pe /= endpe ) then
    call mpp_send ( input, size(input), endpe )
  else
    output = input
    do n=2,npes
      call mpp_recv ( input, size(input), pelist(n) )
      output = output + input
    enddo
  endif

  if ( pe == endpe ) then
    temp = output
    do n=2,npes
      call mpp_send ( temp, size(temp), pelist(n) )
    enddo
  else
    call mpp_recv ( output, size(output), endpe )
  endif

  call mpp_sync_self ( pelist )

  return
end subroutine collect_4d

!#######################################################################

end module zetac_communicate_mod
