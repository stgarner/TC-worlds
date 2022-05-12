subroutine zetac_user_ref_atmos (                                      &
                                              Mgrid, Hmetric, Vmetric, &
                                           uuref, exref, thref, rhref )

use mpp_mod,          only : stdout, mpp_pe, mpp_root_pe
use fms_mod,          only : file_exist, open_namelist_file, close_file
use constants_mod,    only : tfreeze

use zetac_convert_var_mod,       only : get_rh, get_th
use zetac_coriolis_param_mod,    only : coriolis_param
use zetac_horiz_grid_type_mod,   only : horiz_grid_type
use zetac_horiz_metric_type_mod, only : horiz_metric_type
use zetac_phys_con_mod,          only : pref, cpdry, rrat, grav
use zetac_vert_metric_type_mod,  only : vert_metric_type

implicit none

!-----------------------------------------------------------------------
! program to define 2D wind, temperature, pressure and moisture fields
!-----------------------------------------------------------------------

type (horiz_grid_type),                              intent (in)   ::  &
                                                                Mgrid

type (horiz_metric_type),                            intent(in)    ::  &
                                                              Hmetric

type (vert_metric_type),                             intent(in)    ::  &
                                                              Vmetric

real, dimension(Mgrid%jbg:Mgrid%jeg, Mgrid%iz),      intent (out)  ::  &
                                           uuref, exref, thref, rhref

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

real, dimension(Mgrid%iz)                          ::                  &
                                               zeem, exnr, temp, qvap
real, dimension(Mgrid%jbg:Mgrid%jeg, Mgrid%iz)     ::                  &
                                                         taref, ppref
 
!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

real :: edif, gdz
integer :: k, iz

  iz = Mgrid%iz
  zeem = Vmetric%zeem

  call jordan_profiles ( zeem, temp, qvap )

  exnr(1) = 1.0

  do k=2,iz
     gdz = grav*(zeem(k) - zeem(k-1))
     exnr(k) = exnr(k-1)*(1.0 - 0.5*gdz/(cpdry*temp(k-1)))/            &
                         (1.0 + 0.5*gdz/(cpdry*temp(k)))    
  enddo

  edif = 0.5*(exnr(1) + exnr(2)) - 1.0
  exnr = exnr - edif

!-----------------------------------------------------------------------
! generate an axisymmetric, balanced, hurricane-like vortex
!-----------------------------------------------------------------------
 
!-----------------------------------------------------------------------
!  this subroutine generates a balanced axisymmetric, baroclinic vortex. 
!  first we define a reference state at the outer boundary. then we 
!  generate the desired velocity field of the vortex. then we integrate 
!  inwards from the outer boundary to compute the pressure field. 
!  then we update the temperature field to enforce hydrostatic balance.
!  the we repeat these last two steps until the error in the pressure 
!  gradient balance is really small.
!-----------------------------------------------------------------------

  ulat = Hmetric%ulat
  vlat = Hmetric%vlat

  call define_uu ( uuref )
  call init_thref ( thref )

  do i=1,itmax
     call get_exref ( uuref, thref, exref )
     call get_thref ( thref, exref )
     call eval_grad_err ( exref, uuref, thref, errmax )
     if ( i > 1 )then
        if ( abs((errmax-errlast)/errlast) < 0.00001 ) exit
     endif
     errlast = errmax
  enddo

  if ( errlast == errmax ) then
     write ( msg, '('vortex not balanced after',i3,' iterations')' )   &
         itmax
     call error_handler ( msg, WARNING )
  endif

!-----------------------------------------------------------------------
! compute density, potential temperature
!-----------------------------------------------------------------------

  do i=1,nr
     do j=1,nz
        ex(i,j) = (pofrz(i,j)/pref)**(1.0/rkappa)
        rhoofrz  (i,j) = pofrz(i,j)/(rdry*tofrz(i,j))
        thetaofrz(i,j) = tofrz(i,j)/ex(i,j)
     enddo
  enddo

!-----------------------------------------------------------------------
! compute derivative quantities for regularly spaced only
!-----------------------------------------------------------------------

  absvortr = 0.0
  absvortz = 0.0
  dthetadr = 0.0
  dthetadz = 0.0
  
  do i=2,nr-1
    do j=2,nz-1
      dvdz(i,j) = (vofrz(i,j+1) - vofrz(i,j-1))/(2.0*dz)
      dvdr(i,j) = (vofrz(i+1,j) - vofrz(i-1,j))/(2.0*dr)

      absvortz(i,j) = fval + (vofrz(i+1,j) - vofrz(i-1,j))/(2.0*dr)    &
                            + vofrz(i,j)/(i*dr)
      absvortr(i,j) = -dvdz(i,j)

      dthetadr(i,j) = (thetaofrz(i+1,j)-thetaofrz(i-1,j))/(2.0*dr)
      dthetadz(i,j) = (thetaofrz(i,j+1)-thetaofrz(i,j-1))/(2.0*dz)

    enddo
  enddo

!-----------------------------------------------------------------------
! compute symmetric potential vorticity
!-----------------------------------------------------------------------

  do i=1,nr
     do j=1,nz
        qofrz(i,j) = ( absvortr(i,j)*dthetadr(i,j)                     &
                     + absvortz(i,j)*dthetadz(i,j) ) / rhoofrz(i,j)

     enddo
  enddo
  qofrz(1,:) = qofrz(2,:)
  qofrz(nr,:) = qofrz(nr-1,:)
  qofrz(:,1) = qofrz(:,2)
  qofrz(:,nz) = qofrz(:,nz-1)

  do k=1,iz
     taref(:,k) = temp(k)/(1.0 + rrat*qvap(k))
     ppref(:,k) = pref*exnr(k)**rkappa
  enddo

  call get_rh ( ppref, taref, qvref, rhref )
  call get_th ( exref, taref, thref )

  return
end subroutine zetac_user_ref_atmos

!#######################################################################

subroutine eval_grad_err ( pofrz, vofrz, tofrz, z, nr, nz, dr, f, errmax )

real   , dimension (nr,nz), intent (in)  :: pofrz, vofrz, tofrz
real   , dimension (nz),    intent (in)  :: z
integer,                    intent (in)  :: nr, nz
real   ,                    intent (in)  :: dr, f
real   ,                    intent (out) :: errmax

!-----------------------------------------------------------------------
!  local allocations
!-----------------------------------------------------------------------

integer :: i, j
real    :: temp, graderr

  errmax = 0.0

  do i=2,nr-1
     do j=1,nz

        temp = tofrz(i,j)*( 1.0 + 0.608*jordan_q(z(j)) )

!-----------------------------------------------------------------------
! compute error in terms of acceleration
!-----------------------------------------------------------------------

        graderr = ((temp*rdry)/pofrz(i,j))*                            &
                     (pofrz(i+1,j) - pofrz(i-1,j))/(2*dr) -            &
                     (f*vofrz(i,j) + vofrz(i,j)*vofrz(i,j)/((i-1)*dr) )

        errmax  = max (errmax, abs(graderr) )

     enddo
  enddo

  return
end subroutine eval_grad_err

!#######################################################################

subroutine get_exref ( uuref, thref, exref )

real, dimension(Mgrid%jbg:Mgrid%jeg,Mgrid%iz),                         &
                                         intent (in)  :: uuref, thref

real, dimension(Mgrid%jbg:Mgrid%jeg,Mgrid%iz),                         &
                                         intent (out) :: exref

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

integer :: j, k
real    :: fv1, fv2, rad, dr

  iz = Mgrid%iz
  dr = dlat*deg2rad*aearth

  do k=1,iz
     exref(jeg,k) = exnr(k)
     do j=jeg,2,-1
        rad = (ulat(j-1) - rlatmin)*deg2rad*aearth
        fv1 = (fconst*uuref(j-1,k) + uuref(j-1,k)*uuref(j-1,k)/rad)/   &
              (cpdry*thref(j-1,k) )
        rad = (ulat(j) - rlatmin)*deg2rad*aearth
        fv2 = (fconst*uuref(j,k) + uuref(j,k)*uuref(j,k)/rad)/         &
              (cpdry*thref(j,k) )
        exref(j-1,k) = exref(j,k) - 0.5*dr*(fv1 + fv2)
     enddo
  enddo

  return 
end subroutine get_exref

!#######################################################################

subroutine get_thref ( thref, exref )

real   , dimension(nr,nz), intent (in)  :: exref
real   , dimension(nr,nz), intent (out) :: thref

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

integer :: i, j
real    :: dz, dpdz, dz1, dz2, dpdz1, dpdz2, p

!-----------------------------------------------------------------------
! this routine only works for regularly spaced z levels.
!-----------------------------------------------------------------------

  dz = z(2) - z(1)

  do j=1,nz
    do i=1,nr
      p = pofrz(i,j)

      if (j == 1) then

        dpdz = (pofrz(i,3) - pofrz(i,1))/(2.0*dz) -                        &
               dz*(pofrz(i,3) + pofrz(i,1) - 2.0*pofrz(i,2))/(dz*dz) 

      else if (j == nz) then

        dpdz = (pofrz(i,nz) - pofrz(i,nz-2))/(2.0*dz) +                    &
               dz*(pofrz(i,nz) + pofrz(i,nz-2) - 2.0*pofrz(i,nz-1))/(dz*dz) 

      else

        dpdz = (pofrz(i,j+1) - pofrz(i,j-1))/(2.0*dz)
    
      endif    

      tofrz(i,j) = (-p*grav)/(rdry*dpdz*(1 + 0.608*_qvap(k) )
    enddo
  enddo

  return

end subroutine get_thref

!#######################################################################

subroutine define_uu ( uuref )

real, dimension(Mgrid%jbg:Mgrid%jeg,Mgrid%iz), intent (out) :: uuref
 
!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

real, dimension(Mgrid%jbg:Mgrid%jeg) :: uu, vort

integer :: i, j
real    :: rad, dr, ztop, fac, gamma

  gamma =     1.0
  ztop  = 12000.0

!-----------------------------------------------------------------------
! horizontal structure
!-----------------------------------------------------------------------

  call define_vort ( vlat, vort )

  dr = dlat*deg2rad*aearth
  uref(1) = 0.0

  do j=2,jeg
     uu(j) = uu(j-1) + 0.5*dr*vlat(j)*vort(j)
  enddo

  do j=2,jeg
    uu(j) = uu(j)/ulat(j)
  enddo

  do j=jbg,1
     uu(j) = 0.0
  enddo

!-----------------------------------------------------------------------
! vertical structure
!-----------------------------------------------------------------------

  do k=1,iz
     do j=jbg,jeg

        rad = ulat(j)*deg2rad*aearth

!-----------------------------------------------------------------------
!  include next line for vertically aligned profile
!-----------------------------------------------------------------------

        gamma = 0.0

        fac = exp(-((1.7*abs(zee(k))/ztop)**2.0))*                     &
           (1.0 - gamma*exp(-((rad/(15000.0 + 3.0*abs(zee(k)))))**2.7))

        fac = fac*exp( -(rad/350000)**6.0 )

!-----------------------------------------------------------------------
!  if lowest level is below the ground, choose fac so that velocity 
!  field is correct at the surface by linear interpolation
!-----------------------------------------------------------------------

        if (zee(k) < 0.0) then
           fac = 2.0 - fac
        endif

!-----------------------------------------------------------------------
!  include next line to make barotropic vortex
!-----------------------------------------------------------------------

        uuref(j,k) = fac*uu(j)

     enddo
  enddo

  return
end subroutine define_uu

!#######################################################################

subroutine define_vort ( rlat, vort )

real, intent (in)  :: rlat(:)
real               :: vort(:)

real :: rad
integer :: j

!-----------------------------------------------------------------------
! vorticity distribution
!-----------------------------------------------------------------------

  do j=jbg,jeg
     rad = rlat(j)*deg2rad*aearth
     vort(j) = zetamax*exp( -(rad/yscale)**2 )   ! monopolar
  enddo

  return
end subroutine define_vort

!#######################################################################

subroutine jordan_profiles ( zee, temp, qvap )

!-----------------------------------------------------------------------
!  this function returns the temperature as a function of height, based 
!  upon the jordan (1958) mean hurricane season sounding. uses linear 
!  interpolation. will work for altitudes outside of the data range, 
!  uses linear extrapolation
!-----------------------------------------------------------------------

real, intent (in)  :: zee(:)
real, intent (out) :: temp(:), qvap(:)

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

integer :: i, k, kmax
real    :: dz, fac, fac1, t, q

real, dimension(27) :: zdata, tdata, qdata

data zdata /    0.0,   132.0,   583.0,  1054.0,  1547.0,  2063.0,      &
             2609.0,  3182.0,  3792.0,  4442.0,  5138.0,  5888.0,      & 
             6703.0,  7595.0,  8581.0,  9682.0, 10935.0, 12396.0,      & 
            13236.0, 14177.0, 15260.0, 16568.0, 17883.0, 19620.0,      & 
            20743.0, 22139.0, 23971.0 /
    
data tdata /26.3,  26.0,  23.0,  19.8,  17.3,  14.6,  11.8,   8.6,     &
             5.1,   1.4,  -2.5,  -6.9, -11.9, -17.7, -24.8, -33.2,     & 
           -43.3, -55.2, -61.5, -67.6, -72.2, -73.5, -69.8, -63.9,     & 
           -60.6, -57.3, -54.0 /

data qdata / 0.0182, 0.0176, 0.0153, 0.0130, 0.0110, 0.0084, 0.0071,   &
             0.0058, 0.0046, 0.0036, 0.0032, 0.0021, 0.0014, 0.0005,   &
             0.0003, 0.0002, 0.0001, 0.0,    0.0,    0.0,    0.0,      &
             0.0,    0.0,    0.0,    0.0,    0.0,    0.0 /

   kmax = size(zee)
   i = 2

   do k=1,kmax

      do while ( zdata(i) < zee(k) .and. i < 27 )
        i = i + 1
      enddo

      dz = zdata(i) - zdata(i-1)

      fac = (zdata(i) - zee(k))/dz
      fac1 = 1.0 - fac

      t = fac1*tdata(i) + fac*tdata(i-1)
      q = fac1*qdata(i) + fac*qdata(i-1)

      qvap(k) = q*1.1                   ! environmental moistening
      temp(k) = (t + tfreeze)*(1.0 + rrat*qvap(k))

   enddo

   return
end subroutine jordan_profiles

!#######################################################################
   
subroutine error_handler ( message, level ) 
character(len=*), intent(in) :: message
integer,          intent(in) :: level
   
  call error_mesg ( module, message, level)

  return
end subroutine error_handler

!#######################################################################
