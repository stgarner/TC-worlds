module merid_advection_mod ! changed interface

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

public :: merid_advection, merid_advection_end

character(len=128), parameter :: version = '$Id: merid_advection.f90,v 1.1.2.1 2004/08/23 21:11:44 fms Exp $'
character(len=128), parameter :: tagname = '$Name: latest $'

logical :: module_is_initialized = .false.

! buffers for coefficients used by the parabolic scheme
  real, allocatable :: ywts(:,:,:,:), dys(:,:,:)
  integer :: nlons = 0, nlats = 0, nlevs = 0

! for cfl diagnostics with finite volume schemes
  real    :: cflmax = 0.
  integer :: cflerr = 0

contains

!-------------------------------------------------------------------------------

 subroutine merid_advection ( dt, v, dy, r, rdt, mask, scheme, form )

 real, intent(in)                    :: dt
 real, intent(in),  dimension(:,:,:) :: v, dy, r
 real, intent(out), dimension(:,:,:) :: rdt
 real,    intent(in), optional :: mask(:,:,:)
 integer, intent(in), optional :: scheme, form

! see 'vert_advection_mod' for interface documentation

 real, dimension(size(r,1),size(r,2),size(r,3)) :: slp, r_left, r_right
 real, dimension(size(v,1),size(v,2),size(v,3)) :: flux
 real, dimension(0:3,size(r,1),size(r,2),size(r,3)) :: ywt
 real    :: xx, a, b, rm, r6, rst, wt
 real    :: tt, c1, c2
 real    :: small = 1.e-6
 logical :: test_1
 integer :: i, j, k, is, ie, js, je, ks, ke
 integer :: diff_scheme, eqn_form

 real :: cn, rsum, dysum, dtv
 integer :: jj

 if(.not.module_is_initialized) then
   call write_version_number(version, tagname)
   module_is_initialized = .true.
 endif

 ! set default values for optional arguments
   diff_scheme = VAN_LEER_LINEAR
   eqn_form    = FLUX_FORM
   if (present(scheme)) diff_scheme = scheme
   if (present(form))   eqn_form    = form

 ! note: size(r,2)+1 = size(v,2)
   if (size(v,2) /= size(r,2)+1) &
      call error_handler ('horizontal dimension of input arrays inconsistent')

   ks = 1; ke = size(r,3)
   js = 1; je = size(r,2)
   is = 1; ie = size(r,1)

 ! determine fluxes boundaries
 ! most likely v = 0 at these points

   do k = ks, ke
   do i = is, ie
     
   flux(i,js,k)   = v(i,js,k)  *r(i,js,k)
   flux(i,je+1,k) = v(i,je+1,k)*r(i,je,k)
         
   enddo
   enddo
   
   select case (diff_scheme)

   !------ 2nd-order centered scheme assuming variable grid spacing ------
      case (SECOND_CENTERED_WTS)
         do k = ks  , ke
         do j = js+1, je
         do i = is  , ie
            wt = dy(i,j-1,k)/(dy(i,j-1,k)+dy(i,j,k))
            rst = r(i,j-1,k) + wt*(r(i,j,k)-r(i,j-1,k))
            flux(i,j,k) = v(i,j,k)*rst
         enddo
         enddo
         enddo

   !------ 2nd-order centered scheme assumjng uniform grid spacing ------
      case (SECOND_CENTERED)
         do k = ks  , ke
         do j = js+1, je
         do i = is  , ie
            rst = 0.5*(r(i,j,k)+r(i,j-1,k))
            flux(i,j,k) = v(i,j,k)*rst
         enddo
         enddo
         enddo

   !------ 4th-order centered scheme assuming variable grid spacing ------
      case (FOURTH_CENTERED_WTS)
         call compute_weights ( dy, ywt )
         call slope_y ( r, dy, slp, limit=.false., linear=.false. )
         if (present(mask)) then
          ! second order if adiacent to ground
            do k = ks  , ke
            do j = js+2, je-1
            do i = is  , ie
               if (mask(i,j+1,k) > small) then
                  rst = r(i,j-1,k) + ywt(1,i,j,k)*(r(i,j,k)-r(i,j-1,k)) &
                        - ywt(2,i,j,k)*slp(i,j,k) + ywt(3,i,j,k)*slp(i,j-1,k)
               else
                  rst = r(i,j-1,k) + ywt(0,i,j,k)*(r(i,j,k)-r(i,j-1,k))
               endif
            flux(i,j,k) = v(i,j,k)*rst
            enddo
            enddo
            enddo
         else
          ! no mask: always fourth order
            do k = ks  , ke
            do j = js+2, je-1
            do i = is  , ie
               rst = r(i,j-1,k) + ywt(1,i,j,k)*(r(i,j,k)-r(i,j-1,k)) &
                        - ywt(2,i,j,k)*slp(i,j,k) + ywt(3,i,j,k)*slp(i,j-1,k)
               flux(i,j,k) = v(i,j,k)*rst
            enddo
            enddo
            enddo
         endif
         ! second order at top and bottom
         do k = ks, ke
         do i = is, ie
            wt  = dy(i,js,k)/(dy(i,js,k)+dy(i,js+1,k))
            rst = r(i,js,k) + wt*(r(i,js+1,k)-r(i,js,k))
            flux(i,js+1,k) = v(i,js+1,k)*rst
            wt  = dy(i,je-1,k)/(dy(i,je-1,k)+dy(i,je,k))
            rst = r(i,je-1,k) + wt*(r(i,je,k)-r(i,je-1,k))
            flux(i,je,k) = v(i,je,k)*rst
         enddo
         enddo

   !------ 4th-order centered scheme assuming uniform grid spacing ------
      case (FOURTH_CENTERED)
         c1 = 7./12.;  c2 = 1./12.
         if (present(mask)) then
          ! second order if adiacent to ground
            do k = ks  , ke
            do j = js+2, je-1
            do i = is  , ie
               if (mask(i,j+1,k) > small) then
                  rst = c1*(r(i,j,k)+r(i,j-1,k)) - c2*(r(i,j+1,k)+r(i,j-2,k))
               else
                  rst = 0.5*(r(i,j,k)+r(i,j-1,k))
               endif
               flux(i,j,k) = v(i,j,k)*rst
            enddo
            enddo
            enddo
         else
          ! no mask: always fourth order
            do k = ks  , ke
            do j = js+2, je-1
            do i = is  , ie
               rst = c1*(r(i,j,k)+r(i,j-1,k)) - c2*(r(i,j+1,k)+r(i,j-2,k))
               flux(i,j,k) = v(i,j,k)*rst
            enddo
            enddo
            enddo
         endif
         ! second order at top and bottom
         do k = ks, ke
         do i = is, ie
            rst = 0.5*(r(i,js+1,k)+r(i,js,k  ))
            flux(i,js+1,k) = v(i,js+1,k)*rst
            rst = 0.5*(r(i,je,k  )+r(i,je-1,k))
            flux(i,je,k) = v(i,je,k)*rst
         enddo
         enddo

   !------ finite volume scheme using piecewise linear method ------
      case (FINITE_VOLUME_LINEAR)
       ! slope along the z-axis
         call slope_y ( r, dy, slp )
         do k = ks  , ke
         do j = js+1, je
         do i = is  , ie
            if (v(i,j,k) >= 0.) then
               xx = dt*v(i,j,k)/dy(i,j-1,k)
               rst = r(i,j-1,k) + 0.5*slp(i,j-1,k)*(1.-xx)
            else
               xx = -dt*v(i,j,k)/dy(i,j,k)
               rst = r(i,j,k  ) - 0.5*slp(i,j,k  )*(1.-xx)
            endif
            flux(i,j,k) = v(i,j,k)*rst
            if (xx > 1.) cflerr = cflerr+1
            cflmax = max(cflmax,xx)
         enddo
         enddo
         enddo

   !------ finite volume scheme using piecewise parabolic method (PPM) ------
      case (FINITE_VOLUME_PARABOLIC:FINITE_VOLUME_PARABOLIC2)
         call compute_weights ( dy, ywt )
         call slope_y ( r, dy, slp, limit=.true., linear=.false. )
         do k = ks  , ke
         do j = js+2, je-1
         do i = is  , ie
            r_left(i,j,k) = r(i,j-1,k) + ywt(1,i,j,k)*(r(i,j,k)-r(i,j-1,k)) &
                   - ywt(2,i,j,k)*slp(i,j,k) + ywt(3,i,j,k)*slp(i,j-1,k)        
            r_right(i,j-1,k) = r_left(i,j,k)
            ! coming out of this loop, all we need is r_left and r_right
         enddo
         enddo
         enddo

         ! boundary values  ! masks ???????

         do k = ks, ke
         do i = is, ie
           r_left (i,js+1,k) = r(i,js+1,k) - 0.5*slp(i,js+1,k)
           r_right(i,je-1,k) = r(i,je-1,k) + 0.5*slp(i,je-1,k)

           r_right(i,js,k) = r(i,js,k) + 0.5*slp(i,js,k)
           r_left (i,je,k) = r(i,je,k) - 0.5*slp(i,je,k)

           r_left (i,js,k) = r(i,js,k)        ! value not used if v = 0 on boundary
           r_right(i,je,k) = r(i,je,k)        ! value not used if v = 0 on boundary
         enddo
         enddo

     ! monotonicity constraint

       if (diff_scheme == FINITE_VOLUME_PARABOLIC2) then
         ! limiters from Lin (2003), Equation 6 (relaxed constraint)
           do j = js, je
           do k = ks, ke
           do i = is, ie
              r_left (i,j,k) = r(i,j,k) - sign( min(abs(slp(i,j,k)),abs(r_left (i,j,k)-r(i,j,k))), slp(i,j,k) )
              r_right(i,j,k) = r(i,j,k) + sign( min(abs(slp(i,j,k)),abs(r_right(i,j,k)-r(i,j,k))), slp(i,j,k) )
           enddo
           enddo
           enddo
       else
         ! limiters from Colella and Woodward (1984), Equation 1.10
           do j = js, je
           do k = ks, ke
           do i = is, ie
              test_1 = (r_right(i,j,k)-r(i,j,k))*(r(i,j,k)-r_left(i,j,k)) <= 0.0
              if (test_1) then
                 r_left(i,j,k)  = r(i,j,k)
                 r_right(i,j,k) = r(i,j,k)
              endif
              if (j == js .or. j == je) cycle
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
           do j = js+1, je
           do i = is  , ie
              if (v(i,j,k) >= 0.) then
                  cn = dt*v(i,j,k)/dy(i,j-1,k)
                  jj = j-1
                  ! extension for Courant numbers > 1
                  if (cn > 1.) then
                      rsum = 0.
                      dysum = 0.
                      dtv = dt*v(i,j,k)
                      do while (dysum+dy(i,jj,k) < dtv)
                         if (jj == 1) then
                             exit
                         endif
                         dysum = dysum + dy(i,jj,k)
                          rsum =  rsum +  r(i,jj,k)
                         jj = jj-1
                      enddo
                      xx = (dtv-dysum)/dy(i,jj,k)
                  else
                      xx = cn
                  endif
                  rm = r_right(i,jj,k) - r_left(i,jj,k)
                  r6 = 6.0*(r(i,jj,k) - 0.5*(r_right(i,jj,k) + r_left(i,jj,k)))
                  if (jj == js) r6 = 0.
                  rst = r_right(i,jj,k) - 0.5*xx*(rm - (1.0 - tt*xx)*r6)
                  ! extension for Courant numbers > 1
                  if (cn > 1.) rst = (xx*rst + rsum)/cn
              else
                  cn = - dt*v(i,j,k)/dy(i,j,k)
                  jj = j
                  ! extension for Courant numbers > 1
                  if (cn > 1.) then
                      rsum = 0.
                      dysum = 0.
                      dtv = -dt*v(i,j,k)
                      do while (dysum+dy(i,jj,k) < dtv)
                         if (jj == js) then
                             exit
                         endif
                         dysum = dysum + dy(i,jj,k)
                          rsum =  rsum +  r(i,jj,k)
                         jj = jj+1
                      enddo
                      xx = (dtv-dysum)/dy(i,jj,k)
                  else
                      xx = cn
                  endif
                  rm = r_right(i,jj,k) - r_left(i,jj,k)
                  r6 = 6.0*(r(i,jj,k) - 0.5*(r_right(i,jj,k) + r_left(i,jj,k)))
                  if (jj == je) r6 = 0.
                  rst = r_left(i,jj,k) + 0.5*xx*(rm + (1.0 - tt*xx)*r6)
                  ! extension for Courant numbers > 1
                  if (cn > 1.) rst = (xx*rst + rsum)/cn
              endif
              flux(i,j,k) = v(i,j,k)*rst
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
!           rdt (i,j,k) = - (flux(i,j+1,k) - flux (i,j,k)) 
            rdt (i,j,k) = - (flux(i,j+1,k) - flux (i,j,k)) / dy(i,j,k)
       enddo
       enddo
         enddo
      case (ADVECTIVE_FORM)
         do k = ks, ke
         do j = js, je
         do i = is, ie
            rdt (i,j,k) = - (flux(i,j+1,k) - flux (i,j,k) - &
                         r(i,j,k)*(v(i,j+1,k)-v(i,j,k))) / dy(i,j,k)
       enddo
       enddo
         enddo
      case default
        ! ERROR
          call error_handler ('invalid value for optional argument form')
   end select


 end subroutine merid_advection

!-------------------------------------------------------------------------------

 subroutine merid_advection_end

  ! deallocate storage
    if (allocated(ywts)) deallocate(ywts)
    if (allocated(dys))  deallocate(dys)

  ! cfl diagnostics
    call mpp_max (cflmax)
    call mpp_sum (cflerr) ! integer sum
    if (cflmax > 0.) then
        write (stdout(),10) cflmax, cflerr
    endif
 10 format (/,' Merid advection (atmosphere):', &
            /,'     maximum CFL =',f10.6,          &
            /,'     number of CFL errors =',i5,/)
        
 end subroutine merid_advection_end

!-------------------------------------------------------------------------------

 subroutine slope_y ( r, dy, slope, limit, linear )
 real, intent(in),  dimension(:,:,:) :: r, dy
 real, intent(out), dimension(:,:,:) :: slope
 logical, intent(in), optional :: limit, linear

!real    :: grad(2:size(r,1),size(r,2),size(r,3))
 real    :: grad(2:size(r,2))
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
  do i = is, ie

     do j = js+1, je
       grad(j) = (r(i,j,k)-r(i,j-1,k))/(dy(i,j,k)+dy(i,j-1,k))
     enddo
     if (dolinear) then
         do j = js+1, je-1
           slope(i,j,k) = (grad(j+1)+grad(j))*dy(i,j,k)
         enddo
     else
         do j = js+1, je-1
           slope(i,j,k) = (grad(j+1)*(2.*dy(i,j-1,k)+dy(i,j,k)) + &
                           grad(j  )*(2.*dy(i,j+1,k)+dy(i,j,k)))  &
                          *dy(i,j,k)/(dy(i,j-1,k)+dy(i,j,k)+dy(i,j+1,k))
         enddo
     endif
     slope(i,js,k) = 2.*grad(js+1)*dy(i,js,k)
     slope(i,je,k) = 2.*grad(je)  *dy(i,je,k)

   ! apply limiters to slope
     if (limiters) then
        do j = js, je
          if (j >= js+1 .and. j <= je-1) then
            rmin = min(r(i,j-1,k), r(i,j,k), r(i,j+1,k))
            rmax = max(r(i,j-1,k), r(i,j,k), r(i,j+1,k))
            slope(i,j,k) = sign(1.,slope(i,j,k)) *  &
                   min( abs(slope(i,j,k)), 2.*(r(i,j,k)-rmin), 2.*(rmax-r(i,j,k)) )
          else
            slope(i,j,k) = 0.
          endif
        enddo
     endif

  enddo
  enddo

 end subroutine slope_y

!-------------------------------------------------------------------------------

 subroutine compute_weights ( dy, ywt )
 real, intent(in),  dimension(:,:,:)    :: dy
 real, intent(out), dimension(0:,:,:,:) :: ywt
 real    :: denom1, denom2, denom3, denom4, num3, num4, x, y
 integer :: i, j, k, is, ie, js, je, ks, ke
 logical :: redo


   ks = 1; ke = size(dy,3)
   js = 1; je = size(dy,2)
   is = 1; ie = size(dy,1)
 
! check the size of stored coefficients
! need to reallocate if size has changed
   if (nlons /= size(dy,1) .or. nlats /= size(dy,2) .or.  nlevs /= size(dy,3)) then
      if (allocated(ywts)) deallocate(ywts)
      if (allocated(dys))  deallocate(dys)
      nlons = size(dy,1)
      nlats = size(dy,2)
      nlevs = size(dy,3)
      allocate (ywts(0:3,nlons,nlats,nlevs))
      allocate (dys (nlons,nlats,nlevs))
      dys = -1.
   endif
   
! coefficients/weights for computing values at grid box interfaces
! only recompute coefficients for a column when layer depth has changed

   do k = ks, ke
   do i = is, ie

    redo = .false.
    do j = js, je
      if (dy(i,j,k) /= dys(i,j,k)) then
        redo = .true.
        exit
      endif
    enddo

   if (redo) then
     do j = js+2, je-1
       denom1 = 1.0/(dy(i,j-1,k) + dy(i,j,k))
       denom2 = 1.0/(dy(i,j-2,k) + dy(i,j-1,k) + dy(i,j,k) + dy(i,j+1,k))
       denom3 = 1.0/(2*dy(i,j-1,k) +   dy(i,j,k))  
       denom4 = 1.0/(  dy(i,j-1,k) + 2*dy(i,j,k))  
       num3   = dy(i,j-2,k) + dy(i,j-1,k)          
       num4   = dy(i,j,k)   + dy(i,j+1,k)        
       x      = num3*denom3 - num4*denom4        
       y      = 2.0*dy(i,j-1,k)*dy(i,j,k) ! everything up to this point is iust
                                          ! needed to compute x1,x1,x3                      
       ywt(0,i,j,k) = dy(i,j-1,k)*denom1                ! = 1/2 in equally spaced case
       ywt(1,i,j,k) = ywt(0,i,j,k) + x*y*denom1*denom2  ! = 1/2 in equally spaced case
       ywt(2,i,j,k) = dy(i,j-1,k)*num3*denom3*denom2    ! = 1/6 ''
       ywt(3,i,j,k) = dy(i,j,k)*num4*denom4*denom2      ! = 1/6 ''
     enddo
     
     do j = js+2, je-1
       dys (i,j,k)     = dy(i,j,k)
       ywts(0:3,i,j,k) = ywt(0:3,i,j,k)
     enddo
   else

   ! use previously computed coefficients
     do j = js, je
       ywt(0:3,i,j,k) = ywts(0:3,i,j,k)
     enddo
   endif

   enddo
   enddo

 end subroutine compute_weights

!-------------------------------------------------------------------------------

 subroutine error_handler ( message )
 character(len=*), intent(in) :: message

   call error_mesg ('merid_advection', trim(message), FATAL)

 end subroutine error_handler

!-------------------------------------------------------------------------------

end module merid_advection_mod
