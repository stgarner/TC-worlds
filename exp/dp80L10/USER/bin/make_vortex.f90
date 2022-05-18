subroutine zetac_user_perturb (                                        &
                                              Mgrid, Hmetric, Vmetric, &
                                    upert, vpert, epert, tpert, hpert )

use mpp_mod,       only : mpp_pe, mpp_broadcast
use fms_mod,       only : file_exist, open_namelist_file, close_file,  &
                          write_version_number,                        &
                          error_mesg, FATAL, WARNING
use constants_mod, only : tfreeze

use zetac_debug_mod,              only : debug
use zetac_horiz_grid_type_mod,    only : horiz_grid_type
use zetac_horiz_metric_type_mod,  only : horiz_metric_type
use zetac_vert_metric_type_mod,   only : vert_metric_type
use zetac_math_con_mod,           only : pie, deg2rad
use zetac_phys_con_mod,           only : cpdry, rdry, rkappa, rrat,    &
                                         grav, aearth

implicit none

type (horiz_grid_type),   intent(in) ::                         Mgrid

type (horiz_metric_type), intent(in) ::                       Hmetric

type (vert_metric_type),  intent(in) ::                       Vmetric

real, dimension(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed,               &
                Mgrid%kbd:Mgrid%ked),               intent (out)  ::   &
                                    upert, vpert, epert, tpert, hpert

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

integer :: iter, unit, io, length

integer :: i, ibh, ieh, ibg, ieg, ix
integer :: j, jbh, jeh, jbg, jeg, iy
integer :: k, kbd, ked
	   
real :: gdz, thdedz, xhomk
real :: facxy, facz, zee, tvdedz
real :: dlondk, coslat

real, dimension(Mgrid%ibg:Mgrid%ieg) :: ulon
real, dimension(Mgrid%ibg:Mgrid%ieg) :: vlon
real, dimension(Mgrid%jbg:Mgrid%jeg) :: ulat
real, dimension(Mgrid%jbg:Mgrid%jeg) :: vlat

real, dimension(Mgrid%iz) :: zeem

real, dimension(Mgrid%ibg:Mgrid%ieg, Mgrid%jbg:Mgrid%jeg) :: noise

!-----------------------------------------------------------------------
! variables for model fields 
!-----------------------------------------------------------------------

real, dimension(kx) :: zfull
real, dimension(kx) :: zhalf

real, dimension(ix,jx,kx) :: uu
real, dimension(ix,jx,kx) :: vv
real, dimension(ix,jx,kx) :: qr
real, dimension(ix,jx,kx) :: pp
real, dimension(ix,jx,kx) :: th

real, dimension(ix) :: lon
real, dimension(jx) :: lat
real, dimension(ix,jx) :: psurf

real :: dlat, dlon, eval, fval

integer :: i, j, k ,zlev

real :: r, rand

real :: dr = 4000.0
real :: dz =  250.0

real :: dx = 4000.0
real :: dy = 4000.0

real :: cen_lat = 20.05
real :: cen_lon = 55.00

real :: p0 = 101500.0

real :: rb     = 60000.0
real :: zb     =  5000.0
real :: sigmar = 20000.0
real :: sigmaz =  2000.0
real :: pert0  =     0.0

real :: noiseamp   =     1.0
real :: noiserad   = 50000.0
real :: noisewidth = 10000.0

real :: pert, angle
!variable for make_phyd_mod
!variable used by smooth_vort_mod
real, parameter :: zetamax =      .0005
real, parameter :: width   = 60000.

integer, parameter :: itmax=20
real :: errlast, errmax
character(len=*) :: msg

real  ::                    &
  thconst,                  &
  rhconst,                  &
  u0,                       &
  u00,                      &
  htropo,                   &
  dhtropo,                  &
  shtropo,                  &
  bftropo,                  &
  bfstrat,                  &
  rjetsth,                  &
  rjetmid,                  &
  rjetnth

real  ::                    &
  ample=0.0,                &
  xscale=3.0,               &
  yscale=3.0,               &
  zscale=1.0e-3

namelist / reference_nml /  &
  thconst,                  &
  rhconst,                  &
  u0,                       &
  u00,                      &
  htropo,                   &
  dhtropo,                  &
  shtropo,                  &
  bftropo,                  &
  bfstrat,                  &
  rjetsth,                  &
  rjetmid,                  &
  rjetnth 

namelist / perturb_nml /    &
  ample, xscale, yscale, zscale

  do i=1,ix
    do j=1,jx

!-----------------------------------------------------------------------
! mass grid
!-----------------------------------------------------------------------

      xloc = dx*( -ix/2.0 + i )
      yloc = dy*( -jx/2.0 + j )

      do k=1,kx-1
        zloc = zhalf(k)

        call get_data_xyz ( xloc, yloc, zloc, xvel, yvel, pres, temp,  &
                            pofrz, vofrz, tofrz, z, nr, nz, dr)
  
        qr(i,j,k) = jordan_q(zloc)
        pp(i,j,k) = pres
        th(i,j,k) = temp*((100000/pres)**(287.0/1004.0))

!----------------------------------------------------------------------- 
! if desired, add asymmetric bubble or spherical bubble
!-----------------------------------------------------------------------

        r = sqrt((xloc - rb)*(xloc - rb) + yloc*yloc)
                 
        if (xloc .eq. 0.0) then
          if (yloc .ge. 0.0) angle =  pie/2.0
          if (yloc .lt. 0.0) angle = -pie/2.0
        else
          angle = atan(yloc/xloc)
        endif
        if (xloc .lt. 0.0) angle = pie + angle

        pert      = pert0*exp( -((r)/sigmar)**2.0 -((zloc-zb)/sigmaz)**2.0 )
        th(i,j,k) = th(i,j,k) + pert

      enddo

      zloc = 0.0
      call get_data_xyz ( xloc, yloc, zloc, xvel, yvel, pres, temp,    &
                          pofrz, vofrz, tofrz, z, nr, nz, dr )

      psurf(i,j) = pres

!-----------------------------------------------------------------------
! u grid
!-----------------------------------------------------------------------

!      xloc = dx*(-(ix+1)/2.0 + i )
!      yloc = dy*(    -jx/2.0 + j )
      r    = sqrt(xloc*xloc + yloc*yloc)

      do k=1,kx-1
        zloc = zhalf(k)
        call get_data_xyz ( xloc, yloc, zloc, xvel, yvel, pres, temp,  &
                            pofrz, vofrz, tofrz, z, nr, nz, dr)
        uu(i,j,k) = xvel

        call random_number(rand)
	
        uu(i,j,k) = uu(i,j,k) +                                        &
          (rand-0.5)*noiseamp*exp(-((r - noiserad)/noisewidth)**2.0 )* &
                                       exp( - (1.4*zloc/18000)**3.0 )

      enddo

!-----------------------------------------------------------------------
! v grid
!-----------------------------------------------------------------------

!      xloc = dx*(   -(ix)/2.0 + i )
!      yloc = dy*( -(jx+1)/2.0 + j )
      r    = sqrt(xloc*xloc + yloc*yloc)

      do k=1,kx-1
        zloc = zhalf(k)
        call get_data_xyz ( xloc, yloc, zloc, xvel, yvel, pres, temp,  &
                            pofrz, vofrz, tofrz, z, nr, nz, dr )
        vv(i,j,k) = yvel

        call random_number(rand)
    
        vv(i,j,k) = vv(i,j,k) +                                          &
          (rand - 0.5)*noiseamp*exp(-((r - noiserad)/noisewidth)**2.0 )* &
                                         exp( - (1.4*zloc/18000)**3.0 )

      enddo

    enddo
  enddo
