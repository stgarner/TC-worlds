module zetac_horiz_mean_mod

!-----------------------------------------------------------------------
! module information
!-----------------------------------------------------------------------

use mpp_mod,                      only : mpp_pe, mpp_sum, mpp_root_pe, &
                                         mpp_send, mpp_recv
use mpp_domains_mod,              only : mpp_get_layout
use fms_mod,                      only : write_version_number, stdout
use constants_mod,                only : pi

use zetac_horiz_grid_type_mod,    only : horiz_grid_type
use zetac_horiz_metric_type_mod,  only : horiz_metric_type
use zetac_spectral_solver_mod,    only : spectral_solver,              &
                                         spectral_solver_init
use zetac_update_halos_mod,       only : update_halos

!-----------------------------------------------------------------------
! public interfaces
!-----------------------------------------------------------------------

implicit none
private
public horiz_mean, horiz_mean_init

real,    allocatable, dimension(:,:,:) :: sum3
real,    allocatable, dimension(:)     :: wghtx, wghty
real,    allocatable, dimension(:,:,:) :: fx, fy
real,    allocatable, dimension(:,:,:) :: field
real,    allocatable, dimension(:)     :: spng

character(len=*), parameter :: module='zetac_horiz_mean_mod'

real, parameter :: tiny=1.0e-8
integer, parameter :: nvar=4

integer :: ibc, iec, ix, ngridx
integer :: jbc, jec, iy, ngridy
integer :: kbd, ked

integer, allocatable, dimension(:) :: pelist_x
integer, allocatable, dimension(:) :: pelist_y

real :: rad, dx, dy
integer :: pe, endpe_x, endpe_y, npes_x, npes_y

!-----------------------------------------------------------------------
! version control information
!-----------------------------------------------------------------------

character(len=128)  :: version =  '$Id: zetac_horiz_mean.f90,v 1.1.2.8.2.8 2005/07/25 21:13:44 stg Exp $'
character(len=128)  :: tag     =  '$Name:  $'

contains

!#######################################################################

subroutine horiz_mean (                                                &
                                                                Mgrid, &
                                                  facuv, factq, facsp, &
                                           uu_in, vv_in, ta_in, qv_in, &
                       uu_out, vv_out, ta_out, qv_out, ud_out, vd_out )

!-----------------------------------------------------------------------
! gets zonal or global horizontal mean values
!-----------------------------------------------------------------------

type (horiz_grid_type),                               intent (in)  ::  &
                                                                Mgrid

real,                                                  intent(in)  ::  &
                                                  facuv, factq, facsp

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),                 intent (in)  ::  &
                                           uu_in, vv_in, ta_in, qv_in
       
real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),                 intent (out) ::  &
                       uu_out, vv_out, ta_out, qv_out, ud_out, vd_out

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

integer :: i, j, k, nv

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked)                           :: chi

!-----------------------------------------------------------------------
! average in x-y
!-----------------------------------------------------------------------

  do nv=1,nvar

     if (nv == 1) then
        if (facuv == 0.) cycle
        field = uu_in
     else if (nv == 2) then
        if (facuv == 0.) cycle
        field = vv_in
     else if (nv == 3) then
        if (factq == 0.) cycle
        field = ta_in
     else if (nv == 4) then
        if (factq == 0.) cycle
        field = qv_in
     endif

     call get_sum ( Mgrid, field )

      do k=kbd,ked
        do j=jbc,jec
           do i=ibc,iec
              field(i,j,k) = field(i,j,k)/(sum3(i,j,k) + tiny)
           enddo
        enddo    
     enddo

     call update_halos ( Mgrid, field )

     if (nv == 1) then
        uu_out = field
     else if (nv == 2) then
        vv_out = field
     else if (nv == 3) then
        ta_out = field
     else if (nv == 4) then
        qv_out = field
     endif

  enddo  ! loop over variables

!-----------------------------------------------------------------------
! get divergent velocity component
!-----------------------------------------------------------------------

  ud_out = 0.
  vd_out = 0.

  if (facsp == 0.) return

  do k=kbd+1,ked
     do j=jbc,jec
        do i=ibc,iec
           field(i,j,k) = (uu_in(i,j,k) - uu_in(i-1,j,k))/dx           &
                        + (vv_in(i,j,k) - vv_in(i,j-1,k))/dy
        enddo
     enddo
  enddo

!!$  chi = field
!!$  call get_sum ( Mgrid, field )
!!$
!!$  do k=kbd,ked
!!$     do j=jbc,jec
!!$        do i=ibc,iec
!!$           field(i,j,k) = field(i,j,k)/(sum3(i,j,k) + tiny)
!!$        enddo
!!$     enddo
!!$  enddo
!!$
!!$  field = chi - field  ! high-pass
  
  call spectral_solver ( Mgrid, spng, field, chi )
  call update_halos ( Mgrid, chi )

  do k=kbd+1,ked
    if (spng(k) == 0.) cycle
    do j=jbc,jec
       do i=ibc,iec
           ud_out(i,j,k) = (chi(i+1,j,k) - chi(i,j,k))/dx
           vd_out(i,j,k) = (chi(i,j+1,k) - chi(i,j,k))/dy
        enddo
     enddo
  enddo

  ud_out(:,:,kbd) = ud_out(:,:,kbd+1)
  vd_out(:,:,kbd) = vd_out(:,:,kbd+1)

!!$ud_out = uu_in - ud_out
!!$vd_out = vv_in - vd_out

  return
end subroutine horiz_mean

!#######################################################################

subroutine get_sum ( Mgrid, fld )

type (horiz_grid_type),                           intent (in)    ::    &
                                                                Mgrid

real, dimension(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed,               &
                Mgrid%kbd:Mgrid%ked),             intent (inout) ::    &
                                                                  fld
integer :: i, j, k
integer :: i1, i2, j1, j2

  fx = 0.0
  fx(ibc:iec,jbc:jec,:) = fld(ibc:iec,jbc:jec,:)

  call mpp_sum (fx, length=size(fx), pelist=pelist_x)

  do j=jbc,jec

     do k=kbd,ked
        do i=1-ngridx,0
           fx(i,j,k) = fx(i+ix,j,k)
        enddo
        do i=ix+1,ix+ngridx
           fx(i,j,k) = fx(i-ix,j,k)
        enddo
     enddo
     do k=kbd,ked
        do i=ibc,iec
           i1 = i-ngridx
           i2 = i+ngridx
           fld(i,j,k) = sum ( fx(i1:i2,j,k)*wghtx )
        enddo
     enddo

  enddo

  fy = 0.0
  do j=jbc,jec
     do i=ibc,iec
        fy(j,i,:) = fld(i,j,:)
     enddo
  enddo

  call mpp_sum (fy, length=size(fy), pelist=pelist_y)
 
  do i=ibc,iec

     do k=kbd,ked
        do j=1-ngridy,0
           fy(j,i,k) = fy(j+iy,i,k)
        enddo
        do j=iy+1,iy+ngridy
           fy(j,i,k) = fy(j-iy,i,k)
        enddo
     enddo
     do k=kbd,ked
        do j=jbc,jec
           j1 = j-ngridy
           j2 = j+ngridy
           fld(i,j,k) = sum ( fy(j1:j2,i,k)*wghty )
        enddo
     enddo

  enddo

  return
end subroutine get_sum

!#######################################################################

subroutine horiz_mean_init ( Mgrid, Hmetric, rad_in, spfac, spng_in )
 
type (horiz_grid_type),                               intent (in)  ::  &
                                                                Mgrid

type (horiz_metric_type),                             intent (in)  ::  &
                                                              Hmetric
       
real,                                                 intent (in) ::   &
                                                        rad_in, spfac

real, dimension(Mgrid%kbd:Mgrid%ked),                 intent (in) ::   &
                                                              spng_in

real :: arg
integer :: nbuf
integer :: i, j
integer :: layout(2)

!-----------------------------------------------------------------------
! write version number to logfile
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
  nbuf = Mgrid%nbuf

  dx = Hmetric%dx
  dy = Hmetric%dy

  rad = rad_in
  ngridx = int(rad/Hmetric%dx)
  ngridy = int(rad/Hmetric%dy)

  ngridx = max(nbuf, ngridx)
  ngridy = max(nbuf, ngridy)

  write (stdout(), '(/,"zetac_horiz_mean: ngridx ngridy =",2i4)' )   &
                                                         ngridx, ngridy

!
! processor lists
!

  call mpp_get_layout ( Mgrid%Domain, layout )

  npes_x = layout(1)
  npes_y = layout(2)

  allocate ( pelist_x(0:npes_x-1) )
  allocate ( pelist_y(0:npes_y-1) )

  pe = mpp_pe()
  endpe_x = (pe/npes_x)*npes_x
    
  do i=0,npes_x-1
     pelist_x(i) = endpe_x + i
  enddo

  endpe_y = pe - endpe_x
  
  do j=0,npes_y-1
     pelist_y(j) = endpe_y + j*npes_x
  enddo

! Blackman window

  allocate ( wghtx(-ngridx:ngridx) )
  allocate ( wghty(-ngridy:ngridy) )

  wghtx(0) = 1.0
  do i=1,ngridx
     arg = pi*float(i)/ngridx
     wghtx(i) = sqrt( max(0.0,                                       &
                          (0.42 + 0.50*cos(arg) + 0.08*cos(2*arg))) )
!!$     wghtx(i) = max(0.0, 0.42 + 0.50*cos(arg) + 0.08*cos(2*arg))
     wghtx(-i) = wghtx(i)
  enddo
  wghty(0) = 1.0
  do j=1,ngridy
     arg = pi*float(j)/ngridy
     wghty(j) = sqrt( max(0.0,                                       &
                          (0.42 + 0.50*cos(arg) + 0.08*cos(2*arg))) )
!!$     wghty(j) = max(0.0, 0.42 + 0.50*cos(arg) + 0.08*cos(2*arg))
     wghty(-j) = wghty(j)
  enddo

  allocate ( fx(1-ngridx:ix+ngridx,jbc:jec,kbd:ked) )
  allocate ( fy(1-ngridy:iy+ngridy,ibc:iec,kbd:ked) )

  allocate (sum3 (Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed, kbd:ked))
  allocate (field(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed, kbd:ked))

  sum3 = 1.
  call get_sum ( Mgrid, sum3 )

  allocate (spng(kbd:ked))
  spng = spng_in * spfac
  if ( maxval(spng) > 0. ) call spectral_solver_init (Mgrid, Hmetric)

  return
end subroutine horiz_mean_init

!#######################################################################

end module zetac_horiz_mean_mod
