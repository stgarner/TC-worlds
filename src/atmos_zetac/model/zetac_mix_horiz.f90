module zetac_mix_horiz_mod

!-----------------------------------------------------------------------
! module information
!-----------------------------------------------------------------------

use mpp_mod,            only : mpp_pe, mpp_root_pe
use fms_mod,            only : write_version_number,                   &
                               error_mesg, FATAL, WARNING,             &
                               stdlog, stdout
use mpp_domains_mod,    only : mpp_update_domains
use vert_advection_mod, only : ADVECTIVE_FORM, FLUX_FORM

use zetac_extrap_var_mod,        only : extrap_tracer
use zetac_horiz_grid_type_mod,   only : horiz_grid_type
use zetac_horiz_metric_type_mod, only : horiz_metric_type
use zetac_tracer_mod,            only : get_method
use zetac_update_halos_mod,      only : update_halos

!-----------------------------------------------------------------------
! public interfaces
!-----------------------------------------------------------------------

implicit none
private
public mix_horiz, mix_horiz_init

interface mix_horiz
  module procedure mix_horiz_2d
  module procedure mix_horiz_3d
end interface

character(len=*), parameter :: module='zetac_mix_horiz_mod'
character(len=16) :: name_scheme, name_type
integer :: num_form

real :: kh2, kh4                             !  2nd and 4th-order
real :: tqfac

!-----------------------------------------------------------------------
! horizontal diffusion coefficients
!-----------------------------------------------------------------------

real :: c2i0j0, c2i0j1, c2i1j0

character(len=128) :: msg
real, allocatable, dimension(:)     :: dxu
real, allocatable, dimension(:,:)   :: data1
real, allocatable, dimension(:,:,:) :: data

integer :: ibc, iec, ibd, ied, ibc1, iec1
integer :: jbc, jec, jbd, jed, jbc1, jec1
integer :: kbc, kec, kbd, ked
integer :: nbuf

logical :: do_fourth

!-----------------------------------------------------------------------
! version control information
!-----------------------------------------------------------------------

character(len=128)  :: version =  '$Id: zetac_mix_horiz.f90,v 1.1.2.6.2.6 2005/07/25 21:05:40 stg Exp $'
character(len=128)  :: tag     =  '$Name:  $'
     
contains

!#######################################################################

subroutine mix_horiz_3d (                                              &
                                                                Mgrid, &
                                                                  var, &
                                                             mix_tend, &
                                                                index, &
                                                    varref, dp, dpref )

!-----------------------------------------------------------------------
! computes harmonic diffusion
!-----------------------------------------------------------------------

type (horiz_grid_type),                     intent (in)    ::          &
                                                                Mgrid

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),       intent (in)    ::          &
                                                                  var

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),       intent (out)   ::          &
                                                             mix_tend

integer,                                    intent (in)         ::     &
                                                                index

real, optional, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,    &
                Mgrid%kbd:Mgrid%ked),       intent (in)    ::          &
                                                                   dp

real, optional, dimension(Mgrid%jbd:Mgrid%jed, Mgrid%kbd:Mgrid%ked),   &
                                            intent (in)    ::          &
                                                        varref, dpref

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

real    :: factor
integer :: i, j, k
integer :: length
logical :: mix_angmom

  call get_method ( index, 'diffuse_horiz', name_scheme )
  if ( name_scheme(1:4) == 'none' ) return

  do_fourth = ( name_scheme(1:4) == 'four' )
  if ( do_fourth ) then
     factor = -kh4
  else
     factor = kh2
  endif

  call get_method ( index, 'diffuse_type', name_type )

  if ( name_type(1:4) == 'ther' .or. name_type(1:4) == 'trac' ) then
     factor = factor*tqfac
  endif

  if ( factor == 0.0 ) then
     mix_tend = 0.0
     return
  endif

!-----------------------------------------------------------------------
! apply density weighting and subract reference field
!-----------------------------------------------------------------------

  call get_method ( index, 'equation', num=num_form )

  if ( num_form == FLUX_FORM .and. present(dp) ) then
     data = var*dp
     if ( present(varref) .and. present(dpref) ) then
        do i=ibd,ied
           data(i,:,:) = data(i,:,:) - varref(:,:)*dpref(:,:)
        enddo
     endif
  else
     data = var
     if ( present(varref) ) then
        do i=ibd,ied
           data(i,:,:) = data(i,:,:) - varref(:,:)
        enddo
     endif
  endif

!-----------------------------------------------------------------------
! compute second-order harmonic diffusion
!-----------------------------------------------------------------------

  do k=kbc,kec
     do j=jbc1,jec1
        do i=ibc1,iec1
           mix_tend(i,j,k) = ( data(i,j,k)*c2i0j0                      &
            + (data(i,j-1,k) + data(i,j+1,k))*c2i0j1                   &
            + (data(i-1,j,k) + data(i+1,j,k))*c2i1j0 )*factor
        enddo
     enddo
  enddo

!-----------------------------------------------------------------------
! compute fourth-order harmonic diffusion
!-----------------------------------------------------------------------

  if ( do_fourth ) then

     if ( nbuf == 1 ) then
        call extrap_tracer ( Mgrid, mix_tend )
        call update_halos ( Mgrid, mix_tend )
     endif

     data = mix_tend

     do k=kbc,kec
        do j=jbc,jec
           do i=ibc,iec
              mix_tend(i,j,k) = data(i,j,k)*c2i0j0                     &
                             + (data(i,j-1,k) + data(i,j+1,k))*c2i0j1  &
                             + (data(i-1,j,k) + data(i+1,j,k))*c2i1j0
           enddo
        enddo
     enddo

  endif

  if ( num_form == FLUX_FORM .and. present(dp) ) then
     do k=kbc,kec
        do j=jbc,jec
           do i=ibc,iec
              mix_tend(i,j,k) = mix_tend(i,j,k)/dp(i,j,k)
           enddo
        enddo
     enddo
  endif

  if ( mix_angmom ) then
     do j=jbc,jec
        mix_tend(ibc:iec,j,kbc:kec) = mix_tend(ibc:iec,j,kbc:kec)/dxu(j)
     enddo
  endif

  return
end subroutine mix_horiz_3d

!#######################################################################

subroutine mix_horiz_2d (                                              &
                                                                Mgrid, &
                                                                  var, &
                                                             mix_tend, &
                                                                index, &
                                                               varref )

!-----------------------------------------------------------------------
! computes harmonic diffusion
!-----------------------------------------------------------------------

type (horiz_grid_type),                     intent (in)         ::     &
                                                                Mgrid

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed),             &
                                            intent (in)         ::     &
                                                                  var

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed),             &
                                            intent (out)        ::     &
                                                             mix_tend

integer,                                    intent (in)         ::     &
                                                                index

real, optional, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed),   &
                                            intent (in)         ::     &
                                                               varref

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

real    :: factor
integer :: i, j, k
 
  call get_method ( index, 'diffuse_horiz', name_scheme )
  if ( name_scheme(1:4) == 'none' ) return

  do_fourth = ( name_scheme(1:4) == 'four' )
  if ( do_fourth ) then
     factor = -kh4
  else
     factor = kh2
  endif

  call get_method ( index, 'diffuse_type', name_type )

  if ( name_type(1:4) == 'ther' ) then
     factor = factor*tqfac
  endif

  if ( factor == 0.0 ) then
     mix_tend = 0.0
     return
  endif

!-----------------------------------------------------------------------
! subtract reference field
!-----------------------------------------------------------------------

  data1 = var

  if ( present(varref) ) then
     data1 = var - varref
  else
     data1 = var
  endif

!-----------------------------------------------------------------------
! compute second-order harmonic diffusion
!-----------------------------------------------------------------------

  do j=jbc1,jec1
     do i=ibc1,iec1
        mix_tend(i,j) = ( data1(i,j)*c2i0j0                            &
        + (data1(i,j-1) + data1(i,j+1))*c2i0j1                         &
        + (data1(i-1,j) + data1(i+1,j))*c2i1j0 )*factor
     enddo
  enddo

!-----------------------------------------------------------------------
! compute fourth-order harmonic diffusion
!-----------------------------------------------------------------------

  if ( do_fourth ) then

     if ( nbuf == 1 ) then 
        call extrap_tracer ( Mgrid, mix_tend )
        call update_halos ( Mgrid, mix_tend )
     endif
     
     data1 = mix_tend

     do j=jbc,jec
        do i=ibc,iec
           mix_tend(i,j) = data1(i,j)*c2i0j0                           &
                        + (data1(i,j-1) + data1(i,j+1))*c2i0j1         &
                        + (data1(i-1,j) + data1(i+1,j))*c2i1j0
        enddo
     enddo

  endif

  return
end subroutine mix_horiz_2d

!!#######################################################################

subroutine mix_horiz_init (                                            &
                                                       Mgrid, Hmetric, &
                                                           rkh2, rkh4, &
                                                  tqfactor, zonfactor, &
                                                                dtmax )

type (horiz_grid_type),                     intent (in)    ::          &
                                                                Mgrid

type (horiz_metric_type) ,                  intent (in)    ::          &
                                                              Hmetric

real,                                       intent (in)    ::          &
                                                           rkh2, rkh4

real,                                       intent (in)    ::          &
                                           tqfactor, zonfactor, dtmax

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

real    :: gridx, gridy
integer :: io, unit, ierr

  ibc = Mgrid%ibc
  iec = Mgrid%iec
  ibd = Mgrid%ibd
  ied = Mgrid%ied

  jbc = Mgrid%jbc
  jec = Mgrid%jec
  jbd = Mgrid%jbd
  jed = Mgrid%jed

  kbc = Mgrid%kbc
  kec = Mgrid%kec
  kbd = Mgrid%kbd
  ked = Mgrid%ked

!-----------------------------------------------------------------------
! write version number to logfile
!-----------------------------------------------------------------------

  call write_version_number (version, tag)
    
  kh2 = rkh2
  kh4 = rkh4
  tqfac = tqfactor

  if ( Mgrid%iy > 1 ) then
     gridx = Hmetric%dy/zonfactor
     gridy = Hmetric%dy
  else
     gridx = Hmetric%dx
     gridy = Hmetric%dx*zonfactor
  endif

!-----------------------------------------------------------------------
! check for numerical stability
!-----------------------------------------------------------------------

  kh4 = abs(kh4)
  if (kh4  >  0.25*amin1(gridx,gridy)**4/dtmax)  then
     write ( msg, '("kh4=",e15.3," > ",e15.3)' )                       &
                                   kh4, 0.25*amin1(gridx,gridy)**4/dtmax
     call error_handler ( msg )
  endif
  if (kh2 > amin1(gridx,gridy)**2/dtmax)  then
     write ( msg, '("kh2 =",e12.3," >",e12.3)' )                       &
                                       kh2, amin1(gridx,gridy)**2/dtmax
     call error_handler ( msg )
  endif

!-----------------------------------------------------------------------
! coefficients for del^2
!-----------------------------------------------------------------------

  c2i1j0 =  1.0/(gridx*gridx)
  c2i0j1 =  1.0/(gridy*gridy)
  c2i0j0 = -2.0*(c2i1j0 + c2i0j1)

  allocate (data1(ibd:ied,jbd:jed))
  allocate (data (ibd:ied,jbd:jed,Mgrid%kbd:Mgrid%ked))

  data =  0.0 ; data1 = 0.0

  nbuf = Mgrid%nbuf
  if ( nbuf > 1 ) then
     ibc1 = ibc-1
     iec1 = iec+1
     jbc1 = jbc-1
     jec1 = jec+1
  else
     ibc1 = ibc
     iec1 = iec
     jbc1 = jbc
     jec1 = jec
  endif

  allocate (dxu(jbd:jed))
  dxu = Hmetric%dxu(jbd:jed)

  return
end subroutine mix_horiz_init

!#######################################################################

subroutine error_handler ( message ) 
character(len=*), intent(in) :: message
   
  if ( mpp_pe() == mpp_root_pe() ) then
     call error_mesg (module, message, WARNING)
  endif

  return
end subroutine error_handler

!#######################################################################

end module zetac_mix_horiz_mod


