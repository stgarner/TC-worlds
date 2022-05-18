module zetac_spectral_solver_mod

!-----------------------------------------------------------------------
! module information
!-----------------------------------------------------------------------

use mpp_mod, only : mpp_pe, mpp_root_pe, mpp_sum
use fms_mod, only : write_version_number, stdlog, stdout
use fft_mod, only : fft_init, fft_end,                                 &
                    fft_grid_to_fourier,                               &
                    fft_fourier_to_grid
use constants_mod,               only : pi

use zetac_horiz_grid_type_mod,   only : horiz_grid_type
use zetac_horiz_metric_type_mod, only : horiz_metric_type
use zetac_vert_metric_type_mod,  only : vert_metric_type

!-----------------------------------------------------------------------
! public interfaces
!-----------------------------------------------------------------------

implicit none
private
public spectral_solver, spectral_solver_init

character(len=*), parameter :: module='zetac_spectral_solver_mod'

real, parameter :: tiny=1.e-12

!-----------------------------------------------------------------------
! fft arrays
!-----------------------------------------------------------------------

real,             allocatable, dimension(:,:) :: array
real (KIND=4),    allocatable, dimension(:,:) :: fxwave1_r, fxwave1_i
real (KIND=4),    allocatable, dimension(:,:) :: fxwave2_r, fxwave2_i
complex (KIND=4), allocatable, dimension(:,:) :: fxwave1, fxwave2
complex (KIND=4), allocatable, dimension(:,:) :: fxywave1_r, fxywave1_i
complex (KIND=4), allocatable, dimension(:,:) :: fxywave2_r, fxywave2_i

real,      allocatable, dimension(:,:)   :: wght
real,      allocatable, dimension(:)     :: kay, ell

real    :: dx, dy

integer :: ibc, iec, ix, ik
integer :: jbc, jec, iy, il
integer :: kbd, ked

!-----------------------------------------------------------------------
! version control information
!-----------------------------------------------------------------------

character(len=128)  :: version =  '$Id: zetac_spectral_solver.f90,v 1.1.2.4.2.2 2004/07/09 00:47:18 stg Exp $'
character(len=128)  :: tag     =  '$Name: latest $'

contains

!#######################################################################

subroutine spectral_solver ( Mgrid, spng, diverg, velpot )

type (horiz_grid_type),                intent (in)     ::              &
                                                                Mgrid

real, dimension(Mgrid%kbd:Mgrid%ked),  intent (in)     ::              &
                                                                 spng

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),  intent (in)     ::              &
                                                               diverg

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked), intent (out)     ::              &
                                                               velpot

integer :: k

!-----------------------------------------------------------------------
! spectral transform
!-----------------------------------------------------------------------

  do k=kbd+1,ked

     if (spng(k) == 0.) cycle

     array = 0.0
     array(ibc:iec,jbc:jec) = diverg(ibc:iec,jbc:jec,k)

     call mpp_sum ( array, size(array) )

     call spec_transform ( Mgrid, array, velpot(:,:,k) )

  enddo

  return
end subroutine spectral_solver

!#######################################################################

subroutine spec_transform ( Mgrid, fxygrid, chi )

type (horiz_grid_type),                     intent (in)    ::          &
                                                                Mgrid

real, dimension(Mgrid%ix+1, Mgrid%iy),                                 &
                                            intent (inout) ::          &
                                                              fxygrid
      
real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed),             &
                                            intent (out)   ::          &
                                                                  chi

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

integer :: i, j
integer :: k, l

!-----------------------------------------------------------------------
! transform in x from grid to fourier space
!-----------------------------------------------------------------------

  call fft_init (ix)
  fxwave1 = fft_grid_to_fourier ( fxygrid )
  call fft_end

  do j=1,iy
     do k=0,ik
        fxwave1_r(j,k) = REAL (fxwave1(k,j))
        fxwave1_i(j,k) = AIMAG(fxwave1(k,j))
     enddo
  enddo

!-----------------------------------------------------------------------
! transform in y from grid to fourier
!-----------------------------------------------------------------------

  call fft_init (iy)
  fxywave1_r = fft_grid_to_fourier ( fxwave1_r )
  fxywave1_i = fft_grid_to_fourier ( fxwave1_i )

  do l=0,il
     do k=0,ik
        fxywave1_r(l,k) =  wght(k,l)*fxywave1_r(l,k)
        fxywave1_i(l,k) =  wght(k,l)*fxywave1_i(l,k)
     enddo
  enddo

!-----------------------------------------------------------------------
! transform in y from fourier to grid
!-----------------------------------------------------------------------

  fxwave1_r = fft_fourier_to_grid ( fxywave1_r )
  fxwave1_i = fft_fourier_to_grid ( fxywave1_i )
  call fft_end

  do j=1,iy
     do k=0,ik
        fxwave1(k,j) = CMPLX(fxwave1_r(j,k),fxwave1_i(j,k))
     enddo
  enddo
  
!-----------------------------------------------------------------------
! transform in x from fourier to grid
!-----------------------------------------------------------------------

  call fft_init (ix)
  fxygrid = fft_fourier_to_grid ( fxwave1 )
  chi(ibc:iec,jbc:jec) = fxygrid(ibc:iec,jbc:jec)
  call fft_end

  return
end subroutine spec_transform

!#######################################################################

subroutine spectral_solver_init (                                       &
                                                       Mgrid, Hmetric )

type (horiz_grid_type),                     intent (in)    ::          &
                                                                Mgrid
       
type (horiz_metric_type),                   intent (in)    ::          &
                                                              Hmetric

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

integer :: i, j, k, l, l1
real    :: delk, dell

!-----------------------------------------------------------------------
! write version number and namelist to logfile
!-----------------------------------------------------------------------

  call write_version_number (version, tag)

  ibc = Mgrid%ibc
  iec = Mgrid%iec
  jbc = Mgrid%jbc
  jec = Mgrid%jec
  kbd = Mgrid%kbd
  ked = Mgrid%ked

  ix = Mgrid%ix
  iy = Mgrid%iy

  dx = Hmetric%dx
  dy = Hmetric%dy

  allocate (array(ix+1,iy))
  
!-----------------------------------------------------------------------
! define wavenumber arrays and initialize fast fourier transform
!-----------------------------------------------------------------------

  ik = ix/2
  il = iy/2

  write (stdout(),'(/,"zetac_spectral_solver: ik il =",2i4)') ik,il

  allocate (kay(0:ik), ell(0:il))

  delk = 2.0*pi/(ix*dx)
  do k=0,ik
     kay(k) = k*delk
  enddo
  dell = 2.0*pi/(iy*dy)
  do l=0,il
     ell(l) = l*dell
  enddo

  allocate (fxwave1(0:ik,iy))
  allocate (fxwave2(0:ik,iy))
     
  allocate (fxwave1_r(iy+1,0:ik))
  allocate (fxwave1_i(iy+1,0:ik))
  allocate (fxwave2_r(iy+1,0:ik))
  allocate (fxwave2_i(iy+1,0:ik))
     
  allocate (fxywave1_r(0:il,0:ik))
  allocate (fxywave1_i(0:il,0:ik))
  allocate (fxywave2_r(0:il,0:ik))
  allocate (fxywave2_i(0:il,0:ik))

  allocate (wght(0:ik,0:il))

  do l=0,il
     do k=0,ik
        wght(k,l) = -1.0/(kay(k)**2 + ell(l)**2 + tiny)
     enddo
  enddo

  return
end subroutine spectral_solver_init

!#######################################################################

end module zetac_spectral_solver_mod
