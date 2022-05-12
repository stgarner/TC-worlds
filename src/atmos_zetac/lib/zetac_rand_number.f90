module zetac_rand_number_mod

!-----------------------------------------------------------------------
! public interface
!-----------------------------------------------------------------------

implicit none
private
public rand_number

integer, parameter :: ia = 16807
integer, parameter :: ic = 2147483647
integer, parameter :: iq = 127773
integer, parameter :: ir = 2836

integer :: iseed = 123456789

!-----------------------------------------------------------------------
! version control information
!-----------------------------------------------------------------------

character(len=128)  :: version =  '$Id: zetac_rand_number.f90,v 1.1.2.4.2.2 2004/07/09 00:47:18 stg Exp $'
character(len=128)  :: tag     =  '$Name: latest $'

contains

!#######################################################################

subroutine rand_number (r)

!-----------------------------------------------------------------------
!
! uniform random number generator x(n+1) = a*x(n) mod c with
! a = 7**5 and c = 2**(31)-1.
!
! park, steven k. and miller, keith w., "random number generators: 
! good ones are hard to find", communications of the acm, october, 1988.
!
!-----------------------------------------------------------------------

real, intent (out) :: r 

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

integer :: ih, il, it

!-----------------------------------------------------------------------
! initialize seed
!-----------------------------------------------------------------------

  ih = iseed/iq
  il = mod(iseed,iq)
  it = ia*il - ir*ih

  if ( it > 0 ) then
    iseed = it
  else
    iseed = ic + it
  end if

  r = iseed/FLOAT(ic)
  return

end subroutine rand_number

!#######################################################################

end module zetac_rand_number_mod
