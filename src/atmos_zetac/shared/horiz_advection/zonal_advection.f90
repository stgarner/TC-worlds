module zonal_advection_mod ! changed interface

!-------------------------------------------------------------------------------

use fms_mod, only: error_mesg, FATAL, write_version_number, stdout
use mpp_mod, only: mpp_sum, mpp_max, mpp_pe
use vert_advection_mod, only :                                         &
                                                      SECOND_CENTERED, &
                                                      FOURTH_CENTERED, &
                                                 FINITE_VOLUME_LINEAR, &
                                              FINITE_VOLUME_PARABOLIC, &
                                             FINITE_VOLUME_PARABOLIC2, &
                                                  SECOND_CENTERED_WTS, &
                                                  FOURTH_CENTERED_WTS, &
                                VAN_LEER_LINEAR, FINITE_VOLUME_LINEAR, &
                                            FLUX_FORM, ADVECTIVE_FORM

implicit none
private

public :: zonal_advection, zonal_advection_end

character(len=128), parameter :: version = '$Id: zonal_advection.f90,v 1.1.2.1 2004/08/23 21:11:44 fms Exp $'
character(len=128), parameter :: tagname = '$Name: latest $'

logical :: module_is_initialized = .false.

! buffers for coefficients used by the parabolic scheme
  real, allocatable :: xwts(:,:,:,:), dxs(:,:,:)
  integer :: nlons = 0, nlats = 0, nlevs = 0

! for cfl diagnostics with finite volume schemes
  real    :: cflmax = 0.
  integer :: cflerr = 0

contains

!-------------------------------------------------------------------------------

 subroutine zonal_advection ( dt, u, dx, r, rdt, mask, scheme, form )

 real, intent(in)                    :: dt
 real, intent(in),  dimension(:,:,:) :: u, dx, r
 real, intent(out), dimension(:,:,:) :: rdt
 real,    intent(in), optional :: mask(:,:,:)
 integer, intent(in), optional :: scheme, form

! see 'vert_advection_mod' for interface documentation

 real, dimension(size(r,1),size(r,2),size(r,3)) :: slp, r_left, r_right
 real, dimension(size(u,1),size(u,2),size(u,3)) :: flux
 real, dimension(0:3,size(r,1),size(r,2),size(r,3)) :: xwt
 real    :: xx, a, b, rm, r6, rst, wt
 real    :: tt, c1, c2
 real    :: small = 1.e-6
 logical :: test_1
 integer :: i, j, k, is, ie, js, je, ks, ke
 integer :: diff_scheme, eqn_form

 real :: cn, rsum, dxsum, dtu
 integer :: ii

 if(.not.module_is_initialized) then
   call write_version_number(version, tagname)
   module_is_initialized = .true.
 endif

 ! set default values for optional arguments
   diff_scheme = VAN_LEER_LINEAR
   eqn_form    = FLUX_FORM
   if (present(scheme)) diff_scheme = scheme
   if (present(form))   eqn_form    = form

 ! note: size(r,1)+1 = size(u,1)
   if (size(u,1) /= size(r,1)+1) &
      call error_handler ('horizontal dimension of input arrays inconsistent')

   ks = 1; ke = size(r,3)
   js = 1; je = size(r,2)
   is = 1; ie = size(r,1)

 ! determine fluxes boundaries
 ! most likely u = 0 at these points

   do k = ks, ke
   do j = js, je
     
   flux(is,j,k)   = u(is,j,k)  *r(is,j,k)
   flux(ie+1,j,k) = u(ie+1,j,k)*r(ie,j,k)
         
   enddo
   enddo

   select case (diff_scheme)

   !------ 2nd-order centered scheme assuming variable grid spacing ------
      case (SECOND_CENTERED_WTS)
         do k = ks  , ke
         do j = js  , je
         do i = is+1, ie
            wt = dx(i-1,j,k)/(dx(i-1,j,k)+dx(i,j,k))
            rst = r(i-1,j,k) + wt*(r(i,j,k)-r(i-1,j,k))
            flux(i,j,k) = u(i,j,k)*rst
         enddo
         enddo
         enddo

   !------ 2nd-order centered scheme assuming uniform grid spacing ------
      case (SECOND_CENTERED)
         do k = ks  , ke
         do j = js  , je
         do i = is+1, ie
            rst = 0.5*(r(i,j,k)+r(i-1,j,k))
            flux(i,j,k) = u(i,j,k)*rst
         enddo
         enddo
         enddo

   !------ 4th-order centered scheme assuming variable grid spacing ------
      case (FOURTH_CENTERED_WTS)
         call compute_weights ( dx, xwt )
         call slope_x ( r, dx, slp, limit=.false., linear=.false. )
         if (present(mask)) then
          ! second order if adjacent to ground
            do k = ks  , ke
            do j = js  , je
            do i = is+2, ie-1
               if (mask(i+1,j,k) > small) then
                  rst = r(i-1,j,k) + xwt(1,i,j,k)*(r(i,j,k)-r(i-1,j,k)) &
                        - xwt(2,i,j,k)*slp(i,j,k) + xwt(3,i,j,k)*slp(i-1,j,k)
               else
                  rst = r(i-1,j,k) + xwt(0,i,j,k)*(r(i,j,k)-r(i-1,j,k))
               endif
            flux(i,j,k) = u(i,j,k)*rst
            enddo
            enddo
            enddo
         else
          ! no mask: always fourth order
            do k = ks  , ke
            do j = js  , je
            do i = is+2, ie-1
               rst = r(i-1,j,k) + xwt(1,i,j,k)*(r(i,j,k)-r(i-1,j,k)) &
                        - xwt(2,i,j,k)*slp(i,j,k) + xwt(3,i,j,k)*slp(i-1,j,k)
               flux(i,j,k) = u(i,j,k)*rst
            enddo
            enddo
            enddo
         endif
         ! second order at top and bottom
         do k = ks, ke
         do j = js, je
            wt  = dx(is,j,k)/(dx(is,j,k)+dx(is+1,j,k))
            rst = r(is,j,k) + wt*(r(is+1,j,k)-r(is,j,k))
            flux(is+1,j,k) = u(is+1,j,k)*rst
            wt  = dx(ie-1,j,k)/(dx(ie-1,j,k)+dx(ie,j,k))
            rst = r(ie-1,j,k) + wt*(r(ie,j,k)-r(ie-1,j,k))
            flux(ie,j,k) = u(ie,j,k)*rst
         enddo
         enddo

   !------ 4th-order centered scheme assuming uniform grid spacing ------
      case (FOURTH_CENTERED)
         c1 = 7./12.;  c2 = 1./12.
         if (present(mask)) then
          ! second order if adjacent to ground
            do k = ks  , ke
            do j = js  , je
            do i = is+2, ie-1
               if (mask(i+1,j,k) > small) then
                  rst = c1*(r(i,j,k)+r(i-1,j,k)) - c2*(r(i+1,j,k)+r(i-2,j,k))
               else
                  rst = 0.5*(r(i,j,k)+r(i-1,j,k))
               endif
               flux(i,j,k) = u(i,j,k)*rst
            enddo
            enddo
            enddo
         else
          ! no mask: always fourth order
            do k = ks  , ke
            do j = js  , je
            do i = is+2, ie-1
               rst = c1*(r(i,j,k)+r(i-1,j,k)) - c2*(r(i+1,j,k)+r(i-2,j,k))
               flux(i,j,k) = u(i,j,k)*rst
            enddo
            enddo
            enddo
         endif
         ! second order at top and bottom
         do k = ks, ke
         do j = js, je
            rst = 0.5*(r(is+1,j,k)+r(is,j,k  ))
            flux(is+1,j,k) = u(is+1,j,k)*rst
            rst = 0.5*(r(ie,j,k  )+r(ie-1,j,k))
            flux(ie,j,k) = u(ie,j,k)*rst
         enddo
         enddo

   !------ finite volume scheme using piecewise linear method ------
      case (FINITE_VOLUME_LINEAR)
       ! slope along the z-axis
         call slope_x ( r, dx, slp )
         do k = ks  , ke
         do j = js  , je
         do i = is+1, ie
            if (u(i,j,k) >= 0.) then
               xx = dt*u(i,j,k)/dx(i-1,j,k)
               rst = r(i-1,j,k) + 0.5*slp(i-1,j,k)*(1.-xx)
            else
               xx = -dt*u(i,j,k)/dx(i,j,k)
               rst = r(i,j,k  ) - 0.5*slp(i,j,k  )*(1.-xx)
            endif
            flux(i,j,k) = u(i,j,k)*rst
            if (xx > 1.) cflerr = cflerr+1
            cflmax = max(cflmax,xx)
         enddo
         enddo
         enddo

   !------ finite volume scheme using piecewise parabolic method (PPM) ------

      case (FINITE_VOLUME_PARABOLIC:FINITE_VOLUME_PARABOLIC2)

         call compute_weights ( dx, xwt )
         call slope_x ( r, dx, slp, limit=.true., linear=.false. )
         do k = ks  , ke
         do j = js  , je
         do i = is+2, ie-1
            r_left(i,j,k) = r(i-1,j,k) + xwt(1,i,j,k)*(r(i,j,k)-r(i-1,j,k)) &
                   - xwt(2,i,j,k)*slp(i,j,k) + xwt(3,i,j,k)*slp(i-1,j,k)        
            r_right(i-1,j,k) = r_left(i,j,k)
            ! coming out of this loop, all we need is r_left and r_right
         enddo
         enddo
         enddo

         ! boundary values  ! masks ???????

         do k = ks, ke
         do j = js, je
           r_left (is+1,j,k) = r(is+1,j,k) - 0.5*slp(is+1,j,k)
           r_right(ie-1,j,k) = r(ie-1,j,k) + 0.5*slp(ie-1,j,k)

           r_right(is,j,k) = r(is,j,k) + 0.5*slp(is,j,k)
           r_left (ie,j,k) = r(ie,j,k) - 0.5*slp(ie,j,k)

           r_left (is,j,k) = r(is,j,k)        ! value not used if u = 0 on boundary
           r_right(ie,j,k) = r(ie,j,k)        ! value not used if u = 0 on boundary
         enddo
         enddo

     ! monotonicity constraint

       if (diff_scheme == FINITE_VOLUME_PARABOLIC2) then
         ! limiters from Lin (2003), Equation 6 (relaxed constraint)
           do i = is, ie
           do k = ks, ke
           do j = js, je
              r_left (i,j,k) = r(i,j,k) - sign( min(abs(slp(i,j,k)),abs(r_left (i,j,k)-r(i,j,k))), slp(i,j,k) )
              r_right(i,j,k) = r(i,j,k) + sign( min(abs(slp(i,j,k)),abs(r_right(i,j,k)-r(i,j,k))), slp(i,j,k) )
           enddo
           enddo
           enddo
       else
         ! limiters from Colella and Woodward (1984), Equation 1.10
           do k = ks, ke
           do j = js, je
           do i = is, ie
              test_1 = (r_right(i,j,k)-r(i,j,k))*(r(i,j,k)-r_left(i,j,k)) <= 0.0
              if (test_1) then
                 r_left(i,j,k)  = r(i,j,k)
                 r_right(i,j,k) = r(i,j,k)
              endif
              if (i == is .or. i == ie) cycle
              rm = r_right(i,j,k) - r_left(i,j,k)
              a = rm*(r(i,j,k) - 0.5*(r_right(i,j,k) + r_left(i,j,k)))
              b = rm*rm/6.
              if (a >  b) r_left (i,j,k) = 3.0*r(i,j,k) - 2.0*r_right(i,j,k)
              if (a < -b) r_right(i,j,k) = 3.0*r(i,j,k) - 2.0*r_left (i,j,k)
           enddo
           enddo
           enddo
       endif

         ! compute fluxes at interfaces

           tt = 2./3.
           do k = ks  , ke
           do j = js  , je
           do i = is+1, ie
              if (u(i,j,k) >= 0.) then
                  cn = dt*u(i,j,k)/dx(i-1,j,k)
                  ii = i-1
                  ! extension for Courant numbers > 1
                  if (cn > 1.) then
                      rsum = 0.
                      dxsum = 0.
                      dtu = dt*u(i,j,k)
                      do while (dxsum+dx(ii,j,k) < dtu)
                         if (ii == 1) then
                             exit
                         endif
                         dxsum = dxsum + dx(ii,j,k)
                          rsum =  rsum +  r(ii,j,k)
                         ii = ii-1
                      enddo
                      xx = (dtu-dxsum)/dx(ii,j,k)
                  else
                      xx = cn
                  endif
                  rm = r_right(ii,j,k) - r_left(ii,j,k)
                  r6 = 6.0*(r(ii,j,k) - 0.5*(r_right(ii,j,k) + r_left(ii,j,k)))
                  if (ii == is) r6 = 0.
                  rst = r_right(ii,j,k) - 0.5*xx*(rm - (1.0 - tt*xx)*r6)
                  ! extension for Courant numbers > 1
                  if (cn > 1.) rst = (xx*rst + rsum)/cn
              else
                  cn = - dt*u(i,j,k)/dx(i,j,k)
                  ii = i
                  ! extension for Courant numbers > 1
                  if (cn > 1.) then
                      rsum = 0.
                      dxsum = 0.
                      dtu = -dt*u(i,j,k)
                      do while (dxsum+dx(ii,j,k) < dtu)
                         if (ii == is) then
                             exit
                         endif
                         dxsum = dxsum + dx(ii,j,k)
                          rsum =  rsum +  r(ii,j,k)
                         ii = ii+1
                      enddo
                      xx = (dtu-dxsum)/dx(ii,j,k)
                  else
                      xx = cn
                  endif
                  rm = r_right(ii,j,k) - r_left(ii,j,k)
                  r6 = 6.0*(r(ii,j,k) - 0.5*(r_right(ii,j,k) + r_left(ii,j,k)))
                  if (ii == ie) r6 = 0.
                  rst = r_left(ii,j,k) + 0.5*xx*(rm + (1.0 - tt*xx)*r6)
                  ! extension for Courant numbers > 1
                  if (cn > 1.) rst = (xx*rst + rsum)/cn
              endif
              flux(i,j,k) = u(i,j,k)*rst
              if (xx > 1.) cflerr = cflerr+1
              cflmax = max(cflmax,xx)
           enddo
           enddo
           enddo


      case default
        ! ERROR
          call error_handler ('invalid value for optional argument scheme')
   end select


 ! horizontal advective tendency

   select case (eqn_form)
      case (FLUX_FORM)
         do k = ks, ke
         do j = js, je
         do i = is, ie
!           rdt (i,j,k) = - (flux(i+1,j,k) - flux (i,j,k)) 
            rdt (i,j,k) = - (flux(i+1,j,k) - flux (i,j,k)) / dx(i,j,k)
       enddo
       enddo
       enddo

      case (ADVECTIVE_FORM)
         do k = ks, ke
         do i = is, ie
         do j = js, je
            rdt (i,j,k) = - (flux(i+1,j,k) - flux (i,j,k) - &
                         r(i,j,k)*(u(i+1,j,k)-u(i,j,k))) / dx(i,j,k)
       enddo
       enddo
         enddo
      case default
        ! ERROR
          call error_handler ('invalid value for optional argument form')
   end select


 end subroutine zonal_advection

!-------------------------------------------------------------------------------

 subroutine zonal_advection_end

  ! deallocate storage
    if (allocated(xwts)) deallocate(xwts)
    if (allocated(dxs))  deallocate(dxs)

  ! cfl diagnostics
    call mpp_max (cflmax)
    call mpp_sum (cflerr) ! integer sum
    if (cflmax > 0.) then
        write (stdout(),10) cflmax, cflerr
    endif
 10 format (/,' Zonal advection (atmosphere):', &
            /,'     maximum CFL =',f10.6,          &
            /,'     number of CFL errors =',i5,/)
        
 end subroutine zonal_advection_end

!-------------------------------------------------------------------------------

 subroutine slope_x ( r, dx, slope, limit, linear )
 real, intent(in),  dimension(:,:,:) :: r, dx
 real, intent(out), dimension(:,:,:) :: slope
 logical, intent(in), optional :: limit, linear

!real    :: grad(2:size(r,1),size(r,2),size(r,3))
 real    :: grad(2:size(r,1))
 real    :: rmin, rmax
 integer :: i, j, k, is, ie, js, je, ks, ke
 logical :: limiters, dolinear

  limiters = .true.
  if (present(limit))  limiters = limit
  dolinear = .true.
  if (present(linear)) dolinear = linear

   ks = 1; ke = size(r,3)
   js = 1; je = size(r,2)
   is = 1; ie = size(r,1)
  

! compute slope (weighted for unequal levels)

  do k = ks, ke
  do j = js, je

     do i = is+1, ie
       grad(i) = (r(i,j,k)-r(i-1,j,k))/(dx(i,j,k)+dx(i-1,j,k))
     enddo
     if (dolinear) then
         do i = is+1, ie-1
           slope(i,j,k) = (grad(i+1)+grad(i))*dx(i,j,k)
         enddo
     else
         do i = is+1, ie-1
           slope(i,j,k) = (grad(i+1)*(2.*dx(i-1,j,k)+dx(i,j,k)) + &
                           grad(i  )*(2.*dx(i+1,j,k)+dx(i,j,k)))  &
                          *dx(i,j,k)/(dx(i-1,j,k)+dx(i,j,k)+dx(i+1,j,k))
         enddo
     endif
     slope(is,j,k) = 2.*grad(is+1)*dx(is,j,k)
     slope(ie,j,k) = 2.*grad(ie)  *dx(ie,j,k)

   ! apply limiters to slope
     if (limiters) then
        do i = is, ie
          if (i >= is+1 .and. i <= ie-1) then
            rmin = min(r(i-1,j,k), r(i,j,k), r(i+1,j,k))
            rmax = max(r(i-1,j,k), r(i,j,k), r(i+1,j,k))
            slope(i,j,k) = sign(1.,slope(i,j,k)) *  &
                   min( abs(slope(i,j,k)), 2.*(r(i,j,k)-rmin), 2.*(rmax-r(i,j,k)) )
          else
            slope(i,j,k) = 0.
          endif
        enddo
     endif

  enddo
  enddo

 end subroutine slope_x

!-------------------------------------------------------------------------------

 subroutine compute_weights ( dx, xwt )
 real, intent(in),  dimension(:,:,:)    :: dx
 real, intent(out), dimension(0:,:,:,:) :: xwt
 real    :: denom1, denom2, denom3, denom4, num3, num4, x, y
 integer :: i, j, k, is, ie, js, je, ks, ke
 logical :: redo


   ks = 1; ke = size(dx,3)
   js = 1; je = size(dx,2)
   is = 1; ie = size(dx,1)
 
! check the size of stored coefficients
! need to reallocate if size has changed
   if (nlons /= size(dx,1) .or. nlats /= size(dx,2) .or.  nlevs /= size(dx,3)) then
      if (allocated(xwts)) deallocate(xwts)
      if (allocated(dxs))  deallocate(dxs)
      nlons = size(dx,1)
      nlats = size(dx,2)
      nlevs = size(dx,3)
      allocate (xwts(0:3,nlons,nlats,nlevs))
      allocate (dxs (nlons,nlats,nlevs))
      dxs = -1.
   endif
   
! coefficients/weights for computing values at grid box interfaces
! only recompute coefficients for a column when layer depth has changed

   do k = ks, ke
   do j = js, je

    redo = .false.
    do i = is, ie
      if (dx(i,j,k) /= dxs(i,j,k)) then
        redo = .true.
        exit
      endif
    enddo

   if (redo) then
     do i = is+2, ie-1
       denom1 = 1.0/(dx(i-1,j,k) + dx(i,j,k))
       denom2 = 1.0/(dx(i-2,j,k) + dx(i-1,j,k) + dx(i,j,k) + dx(i+1,j,k))
       denom3 = 1.0/(2*dx(i-1,j,k) +   dx(i,j,k))  
       denom4 = 1.0/(  dx(i-1,j,k) + 2*dx(i,j,k))  
       num3   = dx(i-2,j,k) + dx(i-1,j,k)          
       num4   = dx(i,j,k)   + dx(i+1,j,k)        
       x      = num3*denom3 - num4*denom4        
       y      = 2.0*dx(i-1,j,k)*dx(i,j,k) ! everything up to this point is just
                                          ! needed to compute x1,x1,x3                      
       xwt(0,i,j,k) = dx(i-1,j,k)*denom1                ! = 1/2 in equally spaced case
       xwt(1,i,j,k) = xwt(0,i,j,k) + x*y*denom1*denom2  ! = 1/2 in equally spaced case
       xwt(2,i,j,k) = dx(i-1,j,k)*num3*denom3*denom2    ! = 1/6 ''
       xwt(3,i,j,k) = dx(i,j,k)*num4*denom4*denom2      ! = 1/6 ''
     enddo
     
     do i = is+2, ie-1
       dxs (i,j,k)     = dx (i,j,k)
       xwts(0:3,i,j,k) = xwt(0:3,i,j,k)
     enddo
   else

   ! use previously computed coefficients
     do i = is, ie
       xwt(0:3,i,j,k) = xwts(0:3,i,j,k)
     enddo
   endif

   enddo
   enddo

 end subroutine compute_weights

!-------------------------------------------------------------------------------

 subroutine error_handler ( message )
 character(len=*), intent(in) :: message

   call error_mesg ('zonal_advection', trim(message), FATAL)

 end subroutine error_handler

!-------------------------------------------------------------------------------

end module zonal_advection_mod
