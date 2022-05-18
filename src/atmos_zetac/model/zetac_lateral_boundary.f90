module zetac_lateral_boundary_mod

!-----------------------------------------------------------------------
! module information
!-----------------------------------------------------------------------

use mpp_mod, only : mpp_pe, mpp_root_pe, mpp_npes, mpp_sum,            &
                    input_nml_file
use fms_mod, only : file_exist, close_file, open_namelist_file,        &
                    write_version_number, check_nml_error,             &
                    stdlog, stdout

use zetac_horiz_grid_type_mod,   only : horiz_grid_type
use zetac_horiz_metric_mod,      only : CYLINDRICAL
use zetac_horiz_metric_type_mod, only : horiz_metric_type
use zetac_vert_metric_type_mod,  only : vert_metric_type

!-----------------------------------------------------------------------
! public interfaces
!-----------------------------------------------------------------------

implicit none
private
public lateral_boundary, lateral_boundary_init, lateral_signal

character(len=*), parameter :: module='zetac_lateral_boundary_mod'

real, allocatable, dimension(:)     :: dxu, dxv, sums
real, allocatable, dimension(:)     :: cmaxw, cmaxe, cminw, cmine
real, allocatable, dimension(:,:)   :: cpw, cpe
real, allocatable, dimension(:,:)   :: cps, cpn
real, allocatable, dimension(:,:)   :: var1, var2
real, allocatable, dimension(:,:)   :: damp_w, damp_e, damp_s, damp_n
real, allocatable, dimension(:,:)   :: vel_jk, rdz_jk, vel_ik, rdz_ik
real, allocatable, dimension(:,:)   :: dzdyw, dzdye, dzdxs, dzdxn
real, allocatable, dimension(:,:)   :: rdzdyw, rdzdye, rdzdxs, rdzdxn
real, allocatable, dimension(:,:,:) :: varw, vare, vars, varn

real :: dy, delt, dt_save
real :: cmins, cminn, cmaxs, cmaxn

! namelist:
logical :: no_diverge=.false.

real, parameter :: tiny=1.0e-8

integer :: ibc, iec, ibh, ieh, ibu, is, ix
integer :: jbc, jec, jbh, jeh, jbv, js, iy
integer :: kbc, kec

!-----------------------------------------------------------------------
! version control information
!-----------------------------------------------------------------------

character(len=128)  :: version =  '$Id: zetac_lateral_boundary.f90,v 1.1.2.7.2.7 2005/07/05 02:57:53 stg Exp $'
character(len=128)  :: tag     =  '$Name:  $'

contains

!#######################################################################

subroutine lateral_boundary (                                          &
                                                                Mgrid, &
                                                               uu, vv, &
                                                           d_uu, d_vv, &
                                                                uuref )

type (horiz_grid_type),                     intent (in)    ::          &
                                                                Mgrid

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),       intent (in)    ::          &
                                                               uu, vv

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),       intent (inout) ::          &
                                                           d_uu, d_vv

real, dimension(Mgrid%jbd:Mgrid%jed, Mgrid%kbd:Mgrid%ked),             &
                                            intent (in)    ::          &
                                                                uuref

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

real :: ubar, vbar, rdz, sum_vel, sum_rdz
integer :: i, j, k

!-----------------------------------------------------------------------
! semi-implicit boundary tendency using estimated signal speeds
!-----------------------------------------------------------------------

  vel_jk = 0.0
  rdz_jk = 0.0

  if ( ieh > ix ) then
     if ( maxval(cmaxe) > 0.0 ) then
        do k=kbc,kec
           do j=jbc,jec
              vbar = 0.25*(vv(iec,j  ,k) + vv(iec+1,j  ,k)             &
                         + vv(iec,j-1,k) + vv(iec+1,j-1,k))
              js = sign(1.0, vbar)
              vbar = abs(vbar)/dy
              ubar = cpe(j,k)
              d_uu(iec,j,k) = ( d_uu(iec,j,k)                          &
                + vbar*(uu(iec,j-js,k) - uu(iec,j,k)                   &
                      - uuref(j-js,k) + uuref(j,k))                    &
                + ubar*(uu(iec-1,j, k) - uu(iec,j,k)                   &
                      - uuref(j,k) + uuref(j,k)) )/                    &
                 (1.0 + 0.5*delt*(ubar + vbar))
              vel_jk(j,k) = d_uu(iec,j,k)*rdzdye(j,k)
              rdz_jk(j,k) = rdzdye(j,k)
           enddo
        enddo
     else if ( maxval(cmaxe) == 0.0 ) then
        do k=kbc,kec
           do j=jbc,jec
              d_uu(iec,j,k) = 0.0
           enddo
        enddo
     endif
  endif

  if ( ibh < 1 ) then
     if ( maxval(cmaxw) > 0.0 ) then
        do k=kbc,kec
           do j=jbc,jec
              vbar = 0.25*(vv(ibu,j  ,k) + vv(ibu+1,j  ,k)             &
                         + vv(ibu,j-1,k) + vv(ibu+1,j-1,k))
              js = sign(1.0, vbar)
              vbar = abs(vbar)/dy
              ubar = cpw(j,k)
              d_uu(ibu,j,k) = ( d_uu(ibu,j,k)                          &
                + vbar*(uu(ibu,j-js,k) - uu(ibu,j,k)                   &
                      - uuref(j-js,k) + uuref(j,k))                    &
                + ubar*(uu(ibu+1,j, k) - uu(ibu,j,k)                   &
                      - uuref(j,k)  + uuref(j,k)) )/                   &
                 (1.0 + 0.5*delt*(ubar + vbar))
              vel_jk(j,k) = vel_jk(j,k) - d_uu(ibu,j,k)*rdzdyw(j,k)
              rdz_jk(j,k) = rdz_jk(j,k) + rdzdyw(j,k)
           enddo
        enddo
     else if ( maxval(cmaxw) == 0.0 ) then
        do k=kbc,kec
           do j=jbc,jec
              d_uu(ibu,j,k) = 0.0
           enddo
        enddo
     endif
  endif

  vel_ik = 0.0
  rdz_ik = 0.0

  if ( jeh > iy ) then
     if ( cmaxn > 0.0 ) then
        do k=kbc,kec
           do i=ibc,iec
              ubar = 0.25*(uu(i  ,jec,k) + uu(i  ,jec+1,k)             &
                         + uu(i-1,jec,k) + uu(i-1,jec+1,k))
              is = sign(1.0, ubar)
              ubar = abs(ubar)/dxv(jec)
              vbar = cpn(i,k)
              d_vv(i,jec,k) = ( d_vv(i,jec,k)                          &
                + ubar*(vv(i-is,jec,k) - vv(i,jec,k))                  &
                + vbar*(vv(i, jec-1,k) - vv(i,jec,k)) )/               &
                 (1.0 + 0.5*delt*(ubar + vbar))
              vel_ik(i,k) = d_vv(i,jec,k)*rdzdxn(i,k)
              rdz_ik(i,k) = rdzdxn(i,k)
           enddo
        enddo
     else if ( cmaxn == 0.0 ) then
        do k=kbc,kec
           do i=ibc,iec
              d_vv(i,jec,k) = 0.0
           enddo
        enddo
     endif
  endif

  if ( jbh < 1 ) then
     if ( cmaxs > 0.0 ) then
        do k=kbc,kec
           do i=ibc,iec
              ubar = 0.25*(uu(i  ,jbv,k) + uu(i  ,jbv+1,k)             &
                         + uu(i-1,jbv,k) + uu(i-1,jbv+1,k))
              is = sign(1.0, ubar)
              ubar = abs(ubar)/dxv(jbv)
              vbar = cps(i,k)
              d_vv(i,jbv,k) = ( d_vv(i,jbv,k)                          &
                + ubar*(vv(i-is,jbv,k) - vv(i,jbv,k))                  &
                + vbar*(vv(i, jbv+1,k) - vv(i,jbv,k)) )/               &
                 (1.0 + 0.5*delt*(ubar + vbar))
              vel_ik(i,k) = vel_ik(i,k) - d_vv(i,jbv,k)*rdzdxs(i,k)
              rdz_ik(i,k) = rdz_ik(i,k) + rdzdxs(i,k)
           enddo
        enddo
     else if ( cmaxs == 0.0 ) then
        do k=kbc,kec
           do i=ibc,iec
              d_vv(i,jbv,k) = 0.0
           enddo
        enddo
     endif
  endif

!-----------------------------------------------------------------------
! damp away net mass divergence
!-----------------------------------------------------------------------

  if ( no_diverge ) then

     call mpp_sum ( vel_jk, size(vel_jk) )
     call mpp_sum ( rdz_jk, size(rdz_jk) )
     call mpp_sum ( vel_ik, size(vel_ik) )
     call mpp_sum ( rdz_ik, size(rdz_ik) )

     sum_vel = sum( vel_jk ) + sum( vel_ik )
     sum_rdz = sum( rdz_jk ) + sum( rdz_ik )

     if ( sum_rdz > 0.0 ) then

        sum_vel = sum_vel/sum_rdz

        if ( ieh > ix .and. maxval(cmaxe) > 0.0 ) then
           do k=kbc,kec
              do j=jbc,jec
                 d_uu(iec,j,k) = d_uu(iec,j,k) - sum_vel
              enddo
           enddo
        endif
        if ( ibh < 1 .and. maxval(cmaxw) > 0.0 ) then
           do k=kbc,kec
              do j=jbc,jec
                 d_uu(ibu,j,k) = d_uu(ibu,j,k) + sum_vel
              enddo
           enddo
        endif
        if ( jeh > iy .and. cmaxn > 0.0 ) then
           do k=kbc,kec
              do i=ibc,iec
                 d_vv(i,jec,k) = d_vv(i,jec,k) - sum_vel
              enddo
           enddo
        endif
        if ( jbh < 1 .and. cmaxs > 0.0 ) then
           do k=kbc,kec
              do i=ibc,iec
                 d_vv(i,jbv,k) = d_vv(i,jbv,k) + sum_vel
              enddo
           enddo
        endif
     endif

  endif

  return
end subroutine lateral_boundary

!#######################################################################

subroutine lateral_signal (                                            &
                                                                Mgrid, &
                                                          var, varref, &
                                                           nudg, nfac, &
                                               density, qvap, dt_fast )

type (horiz_grid_type),                         intent (in)    ::      &
                                                                Mgrid

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),           intent (in)    ::      &
                                                                  var

real, dimension(Mgrid%jbd:Mgrid%jed, Mgrid%kbd:Mgrid%ked),             &
                                                intent (in)    ::      &
                                                               varref

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),           intent (in)    ::      &
                                                                 nudg

real,                                           intent (in)    ::      &
                                                                 nfac

real, optional, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,    &
                          Mgrid%kbd:Mgrid%ked), intent (in)    ::      &
                                                        density, qvap

real, optional,                                 intent (in)    ::      &
                                                              dt_fast

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

real, dimension(Mgrid%kbd:Mgrid%ked) :: cc
real :: ddx, ddy, ddt
integer :: i, j, k
  
!-----------------------------------------------------------------------
! save fields before start of fast cycle
!-----------------------------------------------------------------------

  if ( present(dt_fast) ) then

     if ( ibh < 1.and. maxval(cmaxw) > 0.0 ) then
        do k=kbc,kec
           do j=jbc,jec
              varw(1,j,k) = var(1,j,k) - varref(j,k)
              varw(2,j,k) = var(2,j,k) - varref(j,k)
              varw(3,j,k) = var(3,j,k) - varref(j,k)
              damp_w(j,k) = nfac*nudg(2,j,k)*varw(2,j,k)
              rdzdyw(j,k) = density(ibc,j,k)*dzdyw(j,k)*              &
                        (1.0 - qvap(ibc,j,k))
           enddo
        enddo
     endif
     if ( ieh > ix .and. maxval(cmaxe) > 0.0 ) then
        do k=kbc,kec
           do j=jbc,jec
              vare(1,j,k) = var(iec  ,j,k) - varref(j,k)
              vare(2,j,k) = var(iec-1,j,k) - varref(j,k)
              vare(3,j,k) = var(iec-2,j,k) - varref(j,k)
              damp_e(j,k) = nfac*nudg(iec-1,j,k)*vare(2,j,k)
              rdzdye(j,k) = density(iec,j,k)*dzdye(j,k)*               &
                        (1.0 - qvap(iec,j,k))
           enddo
        enddo
     endif
     if ( jbh < 1 .and. cmaxs > 0.0 ) then
        do k=kbc,kec
           do i=ibc,iec
              vars(1,i,k) = var(i,1,k) - varref(1,k)
              vars(2,i,k) = var(i,2,k) - varref(2,k)
              vars(3,i,k) = var(i,3,k) - varref(3,k)
              damp_s(i,k) = nfac*nudg(i,2,k)*vars(2,i,k)
              rdzdxs(i,k) = density(i,jbc,k)*dzdxs(i,k)*              &
                        (1.0 - qvap(i,jbc,k))
           enddo
        enddo
     endif
     if ( jeh > iy .and. cmaxn > 0.0 ) then
        do k=kbc,kec
           do i=ibc,iec
              varn(1,i,k) = var(i,jec  ,k) - varref(jec  ,k)
              varn(2,i,k) = var(i,jec-1,k) - varref(jec-1,k)
              varn(3,i,k) = var(i,jec-2,k) - varref(jec-2,k)
              damp_n(i,k) = nfac*nudg(i,jec-1,k)*varn(2,i,k)
              rdzdxn(i,k) = density(i,jec,k)*dzdxn(i,k)*               &
                        (1.0 - qvap(i,jec,k))
           enddo
        enddo
     endif

     dt_save = dt_fast
     delt = 0.5*dt_fast

     return

  endif

!----------------------------------------------------------------------- 
! set fast time step
!-----------------------------------------------------------------------

  cpw = 0.0
  cpe = 0.0
  cpn = 0.0
  cps = 0.0

!----------------------------------------------------------------------- 
! zonal signal speeds at western boundary
!-----------------------------------------------------------------------

  if ( ibh < 1 .and. maxval(cmaxw) > 0.0 ) then

     do j=jbc,jec

        do i=1,3
           do k=kbc,kec
              var1(i,k) = varw(i,j,k)
              var2(i,k) = var(i,j,k) - varref(j,k)
           enddo
        enddo

        do k=kbc,kec
           ddx = 0.25*(var1(3,k) + var2(3,k) - (var1(1,k) + var2(1,k)))
           ddt = (var2(2,k) - var1(2,k)) + delt*damp_w(j,k)
           cc(k) = min(cmaxw(j), max(cminw(j), ddt/(ddx + tiny) ))
        enddo

        do k=kbc+1,kec-1
           cpw(j,k) = (cc(k-1) + cc(k) + cc(k+1))/(3.0*delt)
        enddo

        cpw(j,kbc) = (2.0*cc(kbc) + cc(kbc+1))/(3.0*delt)
        cpw(j,kec) = (2.0*cc(kec) + cc(kec-1))/(3.0*delt)

        do i=1,3
           do k=kbc,kec
              varw(i,j,k) = var2(i,k)
           enddo
        enddo

     enddo
      
  endif

!----------------------------------------------------------------------- 
! zonal signal speeds at eastern boundary
!-----------------------------------------------------------------------

  if ( ieh > ix .and. maxval(cmaxe) > 0.0 ) then

     do j=jbc,jec

        do i=1,3
           do k=kbc,kec
              var1(i,k) = vare(i,j,k)
              var2(i,k) = var(iec+1-i,j,k) - varref(j,k)
           enddo
        enddo

        do k=kbc,kec
           ddx = 0.25*(var1(3,k) + var2(3,k) - (var1(1,k) + var2(1,k)))
           ddt = (var2(2,k) - var1(2,k)) + delt*damp_e(j,k)
           cc(k) = min(cmaxe(j), max(cmine(j), ddt/(ddx + tiny) ))
        enddo

        do k=kbc+1,kec-1
           cpe(j,k) = (cc(k-1) + cc(k) + cc(k+1))/(3.0*delt)
        enddo

        cpe(j,kbc) = (2.0*cc(kbc) + cc(kbc+1))/(3.0*delt)
        cpe(j,kec) = (2.0*cc(kec) + cc(kec-1))/(3.0*delt)

        do i=1,3
           do k=kbc,kec
              vare(i,j,k) = var2(i,k)
           enddo
        enddo

     enddo

  endif

!-----------------------------------------------------------------------
! meridional signal speeds at southern boundary
!-----------------------------------------------------------------------

  if ( jbh < 1 .and. cmaxs > 0.0 ) then

     do i=ibc,iec

        do j=1,3
           do k=kbc,kec
              var1(j,k) = vars(j,i,k)
              var2(j,k) = var(i,j,k) - varref(j,k)
           enddo
        enddo

        do k=kbc,kec
           ddy = 0.25*(var1(3,k) + var2(3,k) - (var1(1,k) + var2(1,k)))
           ddt = (var2(2,k) - var1(2,k)) + delt*damp_s(i,k)
           cc(k) = min(cmaxs, max(cmins, ddt/(ddy + tiny) ))
        enddo

        do k=kbc+1,kec-1
           cps(i,k) = (cc(k-1) + cc(k) + cc(k+1))/(3.0*delt)
        enddo

        cps(i,kbc) = (2.0*cc(kbc) + cc(kbc+1))/(3.0*delt)
        cps(i,kec) = (2.0*cc(kec) + cc(kec-1))/(3.0*delt)

        do j=1,3
           do k=kbc,kec
              vars(j,i,k) = var2(j,k)
           enddo
        enddo

     enddo

  endif

!-----------------------------------------------------------------------
! meridional signal speeds at northern boundary
!-----------------------------------------------------------------------

  if ( jeh > iy .and. cmaxn > 0.0 ) then

     do i=ibc,iec

        do j=1,3
           do k=kbc,kec
              var1(j,k) = varn(j,i,k)
              var2(j,k) = var(i,jec+1-j,k) - varref(jec+1-j,k)
           enddo
        enddo

        do k=kbc,kec
           ddy = 0.25*(var1(3,k) + var2(3,k) - (var1(1,k) + var2(1,k)))
           ddt = (var2(2,k) - var1(2,k)) + delt*damp_n(i,k)
           cc(k) = min(cmaxn, max(cminn, ddt/(ddy + tiny) ))
        enddo

        do k=kbc+1,kec-1
           cpn(i,k) = (cc(k-1) + cc(k) + cc(k+1))/(3.0*delt)
        enddo

        cpn(i,kbc) = (2.0*cc(kbc) + cc(kbc+1))/(3.0*delt)
        cpn(i,kec) = (2.0*cc(kec) + cc(kec-1))/(3.0*delt)

        do j=1,3
           do k=kbc,kec
              varn(j,i,k) = var2(j,k)
           enddo
        enddo

     enddo

  endif

  delt = dt_save

  return
end subroutine lateral_signal

!#######################################################################

subroutine lateral_boundary_init ( Mgrid, Hmetric, Vmetric, dtmax )
 
type (horiz_grid_type),                     intent (in)    ::          &
                                                                Mgrid

type (horiz_metric_type),                   intent (in)    ::          &
                                                              Hmetric

type (vert_metric_type),                    intent (in)    ::          &
                                                              Vmetric

real,                                       intent (in)    ::          &
                                                                dtmax

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

integer :: io, unit, ierr
integer :: i, j, k
    
real ::                           cmin_w, cmin_e, cmin_s, cmin_n,      &
                                  cmax_w, cmax_e, cmax_s, cmax_n

namelist /lateral_boundary_nml/   cmin_w, cmin_e, cmin_s, cmin_n,      &
                                  cmax_w, cmax_e, cmax_s, cmax_n,      &
                                  no_diverge

  ix  = Mgrid%ix
  iy  = Mgrid%iy

  ibc = Mgrid%ibc
  iec = Mgrid%iec
  ibh = Mgrid%ibh
  ieh = Mgrid%ieh
  ibu = Mgrid%ibu

  jbc = Mgrid%jbc
  jec = Mgrid%jec
  jbh = Mgrid%jbh
  jeh = Mgrid%jeh
  jbv = Mgrid%jbv

  kbc = Mgrid%kbc
  kec = Mgrid%kec

!-----------------------------------------------------------------------
! read namelist
!-----------------------------------------------------------------------

  read (input_nml_file, nml=lateral_boundary_nml, iostat=io)
  ierr = check_nml_error(io,'lateral_boundary_nml')
 
!-----------------------------------------------------------------------
! write version number and namelist to logfile
!-----------------------------------------------------------------------

  call write_version_number (version, tag)
  if (mpp_pe() == mpp_root_pe())                                       &
                             write (stdlog(), nml=lateral_boundary_nml)
 
  allocate (cpw(Mgrid%jbd:Mgrid%jed,Mgrid%kbd:Mgrid%ked))
  allocate (cpe(Mgrid%jbd:Mgrid%jed,Mgrid%kbd:Mgrid%ked))
  allocate (cps(Mgrid%ibd:Mgrid%ied,Mgrid%kbd:Mgrid%ked))
  allocate (cpn(Mgrid%ibd:Mgrid%ied,Mgrid%kbd:Mgrid%ked))

  allocate (damp_w(Mgrid%jbd:Mgrid%jed,Mgrid%kbd:Mgrid%ked))
  allocate (damp_e(Mgrid%jbd:Mgrid%jed,Mgrid%kbd:Mgrid%ked))
  allocate (damp_s(Mgrid%ibd:Mgrid%ied,Mgrid%kbd:Mgrid%ked))
  allocate (damp_n(Mgrid%ibd:Mgrid%ied,Mgrid%kbd:Mgrid%ked))
  
  allocate (varw(3,Mgrid%jbd:Mgrid%jed,Mgrid%kbd:Mgrid%ked))
  allocate (vare(3,Mgrid%jbd:Mgrid%jed,Mgrid%kbd:Mgrid%ked))
  allocate (vars(3,Mgrid%ibd:Mgrid%ied,Mgrid%kbd:Mgrid%ked))
  allocate (varn(3,Mgrid%ibd:Mgrid%ied,Mgrid%kbd:Mgrid%ked))
  allocate (var1(3,Mgrid%kbd:Mgrid%ked))
  allocate (var2(3,Mgrid%kbd:Mgrid%ked))

  allocate (dxu(jbc:jec))
  allocate (dxv(jbv:jec))
  
  dxu = Hmetric%dxu(jbc:jec)
  dxv = Hmetric%dxv(jbv:jec)
  dy  = Hmetric%dy

  allocate (cmaxw(jbc:jec))
  allocate (cmaxe(jbc:jec))
  allocate (cminw(jbc:jec))
  allocate (cmine(jbc:jec))

  cmaxw = min(1.0,   cmax_w/dxu*dtmax)
  cmaxe = min(1.0,   cmax_e/dxu*dtmax)
  cminw = min(cmaxw, cmin_w/dxu*dtmax)
  cmine = min(cmaxe, cmin_e/dxu*dtmax)

  cmaxs = min(1.0,   cmax_s/dy*dtmax)
  cmaxn = min(1.0,   cmax_n/dy*dtmax)
  cmins = min(cmaxs, cmin_s/dy*dtmax)
  cminn = min(cmaxn, cmin_n/dy*dtmax)

  if ( ix == 1 ) then
     cmaxw = -1.0
     cmaxe = -1.0
  endif
  if ( iy == 1 ) then
     cmaxs = -1.0
     cmaxn = -1.0
  endif
  if ( Hmetric%geometry == CYLINDRICAL ) then
     cmaxs = 0.0
  endif

  allocate ( dzdyw(Mgrid%jbd:Mgrid%jed,Mgrid%kbd:Mgrid%ked))
  allocate ( dzdye(Mgrid%jbd:Mgrid%jed,Mgrid%kbd:Mgrid%ked))
  allocate ( dzdxs(Mgrid%ibd:Mgrid%ied,Mgrid%kbd:Mgrid%ked))
  allocate ( dzdxn(Mgrid%ibd:Mgrid%ied,Mgrid%kbd:Mgrid%ked))
  allocate (rdzdyw(Mgrid%jbd:Mgrid%jed,Mgrid%kbd:Mgrid%ked))
  allocate (rdzdye(Mgrid%jbd:Mgrid%jed,Mgrid%kbd:Mgrid%ked))
  allocate (rdzdxs(Mgrid%ibd:Mgrid%ied,Mgrid%kbd:Mgrid%ked))
  allocate (rdzdxn(Mgrid%ibd:Mgrid%ied,Mgrid%kbd:Mgrid%ked))

  do k=kbc,kec
     do j=jbc,jec
        dzdyw(j,k) = Vmetric%dzetam(k)*dy
        dzdye(j,k) = Vmetric%dzetam(k)*dy
     enddo
     do i=ibc,iec
        dzdxs(i,k) = Vmetric%dzetam(k)*dxv(jbv)
        dzdxn(i,k) = Vmetric%dzetam(k)*dxv(jec)
     enddo
  enddo

  allocate (vel_ik(ix,kbc:kec))
  allocate (rdz_ik(ix,kbc:kec))
  allocate (vel_jk(iy,kbc:kec))
  allocate (rdz_jk(iy,kbc:kec))

  return
end subroutine lateral_boundary_init

!#######################################################################

end module zetac_lateral_boundary_mod
