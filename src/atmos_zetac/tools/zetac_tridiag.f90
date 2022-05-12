module zetac_tridiag_mod

!-----------------------------------------------------------------------
! module information
!-----------------------------------------------------------------------

use fms_mod, only : write_version_number

!-----------------------------------------------------------------------
! public interfaces
!-----------------------------------------------------------------------
implicit none
private
public solve_tridiag

character(len=*), parameter :: module='zetac_tridiag_mod'
logical :: do_init=.true.

!-----------------------------------------------------------------------
! version control information
!-----------------------------------------------------------------------

character(len=128)  :: version =  '$Id: zetac_tridiag.f90,v 1.1.2.3 2004/05/21 15:46:03 ck Exp $'
character(len=128)  :: tag     =  '$Name: latest $'

contains

!#######################################################################

subroutine solve_tridiag ( bet, gam, a, b, c, rhs, ans )

!-----------------------------------------------------------------------
! solves the m tridiagonal matrix equations au=r, where a is an nxn 
! matrix with lower diagonal a(n), main diagonal b(n), and upper 
! diagonal c(n), solution vector ans(id,n), and right hand side rhs(id,n).
! elements a(*,1) and c(*,n) are not referenced.
!-----------------------------------------------------------------------

real, intent (inout), dimension(:,:,:) :: bet, gam
real, intent (in),    dimension(:,:,:) :: a
real, intent (in),    dimension(:,:,:), optional :: b, c, rhs
real, intent (out),   dimension(:,:,:), optional :: ans

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

integer :: i, id
integer :: j, jd
integer :: k, kd

!-----------------------------------------------------------------------
! initialize module
!-----------------------------------------------------------------------

  if (do_init) then
     call tridiag_init
  endif

  id = size(a,1)
  jd = size(a,2)
  kd = size(a,3)

  if ( present(b) ) then

!-----------------------------------------------------------------------
! perform decomposition
!-----------------------------------------------------------------------

     do j=1,jd
        do i=1,id
           bet(i,j,1) = 1.0/b(i,j,1)     ! b is positive definite
        enddo
     enddo

     do k=2,kd
        do j=1,jd
           do i=1,id
              gam(i,j,k) = c(i,j,k-1)*bet(i,j,k-1)
              bet(i,j,k) = 1.0/(b(i,j,k) - a(i,j,k)*gam(i,j,k))
           enddo
        enddo
     enddo

  endif

  if ( present(rhs) ) then

!-----------------------------------------------------------------------
! perform forward substitution
!-----------------------------------------------------------------------

     do j=1,jd
        do i=1,id
           ans(i,j,1) = rhs(i,j,1)*bet(i,j,1)
        enddo
     enddo

     do k=2,kd
        do j=1,jd
           do i=1,id
              ans(i,j,k) = (rhs(i,j,k) - a(i,j,k)*ans(i,j,k-1))*       &
                                                             bet(i,j,k)
           enddo
        enddo
     enddo

!-----------------------------------------------------------------------
! perform back substitution
!-----------------------------------------------------------------------

     do k=kd-1,1,-1
        do j=1,jd
           do i=1,id
              ans(i,j,k) = ans(i,j,k) - gam(i,j,k+1)*ans(i,j,k+1)
           enddo
        enddo
     enddo

  endif

  return
end subroutine solve_tridiag

!#######################################################################

subroutine tridiag_init

!-----------------------------------------------------------------------
! write version number to logfile
!-----------------------------------------------------------------------

  call write_version_number (version, tag)

  do_init = .false.
  
  return
end subroutine tridiag_init

!#######################################################################

end module zetac_tridiag_mod
