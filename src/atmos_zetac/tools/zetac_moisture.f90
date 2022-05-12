module zetac_moisture_mod

!-----------------------------------------------------------------------
! module information
!-----------------------------------------------------------------------

use mpp_mod,            only : input_nml_file
use fms_mod,            only : write_version_number, check_nml_error
use constants_mod,      only : rdgas, rvgas, hlv, tfreeze

use zetac_math_con_mod, only : ln10

!-----------------------------------------------------------------------
! public interfaces
!-----------------------------------------------------------------------

implicit none
private
public get_qsat, get_esat, get_dlnesatdt, moisture_init

interface get_esat
  module procedure get_esat_3d
  module procedure get_esat_2d
  module procedure get_esat_0d
end interface

interface get_qsat
  module procedure get_qsat_3d
  module procedure get_qsat_2d
  module procedure get_qsat_0d
end interface

character(len=*), parameter :: module='zetac_moisture_mod'
logical :: do_init=.true.
!logical, parameter :: do_CC=.false.

real :: C1, C2, C3, C4, C5, C6, C7, C8, D1, D2
real :: E1, E2, F1, F2, F3, F4, F5

  DATA C1/-9.09718/, C2/-3.56654/, C3/ 0.876793 /, C4/2.78583503/
  DATA C5/-7.90298/, C6/ 5.02808/, C7/-1.3816E-7/, C8/0.0081328/
  DATA D1/11.344  /, D2/-3.49149/
  DATA F5/5.0057149/

real, parameter :: tboil=tfreeze+100.0
real, parameter :: atet=611.0, btet=17.27
real :: logtf, logtb, dtemp, T1

real, parameter :: eps = rdgas/rvgas
real, parameter :: d378 = 1.0-eps
real, parameter :: T0 = 250., Trange = 60.
real, parameter :: eref = 4.55522
integer, parameter :: num = 120

real, allocatable :: etable(:)

!-----------------------------------------------------------------------
! namelist parameters
!-----------------------------------------------------------------------

logical :: do_kessler=.true.
logical :: do_ras=.false.
logical :: do_CC=.false.
real :: qvmin=0.

!-----------------------------------------------------------------------
! version control information
!-----------------------------------------------------------------------

character(len=128)  :: version =  '$Id: zetac_moisture.f90,v 1.1.2.5.2.1 2004/11/04 22:46:16 stg Exp $'
character(len=128)  :: tag     =  '$Name:  $'

contains

!#######################################################################

subroutine get_qsat_3d ( pres, tabs, qsat, wice )

!-----------------------------------------------------------------------
! computes saturation specific humidity
!-----------------------------------------------------------------------

real, dimension(:,:,:),           intent (in)    ::    pres, tabs
real, dimension(:,:,:), optional, intent (in)    ::    wice
real, dimension(:,:,:),           intent (out)   ::    qsat

real, dimension(size(tabs,1),size(tabs,2),size(tabs,3)) :: esat

  if (do_init) then
     call moisture_init
  endif

  call get_esat_3d ( tabs, esat, wice )
  
!  qsat = eps*esat/(pres - d378*esat)
  qsat = eps*esat/(pres - esat) ! mixing ratio

  return
end subroutine get_qsat_3d

!#######################################################################

subroutine get_qsat_2d ( pres, tabs, qsat, wice )

!-----------------------------------------------------------------------
! computes saturation specific humidity
!-----------------------------------------------------------------------

real, dimension(:,:),           intent (in)    ::   pres, tabs
real, dimension(:,:), optional, intent (in)    ::   wice
real, dimension(:,:),           intent (out)   ::   qsat

real, dimension(size(tabs,1),size(tabs,2)) :: esat
integer :: j, k

  if (do_init) then
     call moisture_init
  endif

  call get_esat_2d ( tabs, esat, wice )

!  qsat = eps*esat/(pres - d378*esat)
  qsat = eps*esat/(pres - esat) ! mixing ratio

  return
end subroutine get_qsat_2d

!#######################################################################

subroutine get_qsat_0d ( pres, tabs, qsat )

!-----------------------------------------------------------------------
! computes saturation specific humidity
!-----------------------------------------------------------------------

real, intent (in)  :: pres, tabs
real, intent (out) :: qsat

real :: esat

  if (do_init) then
     call moisture_init
  endif

  call get_esat_0d ( tabs, esat )

  qsat = eps*esat/(pres - esat) ! mixing ratio

  return
end subroutine get_qsat_0d

!#######################################################################

subroutine get_esat_3d ( tabs, esat, wice )

real, dimension(:,:,:), intent (in)  :: tabs
real, dimension(:,:,:), optional, intent (in)  :: wice
real, dimension(:,:,:), intent (out) :: esat

real, dimension(size(tabs,1),size(tabs,2),size(tabs,3)) :: logtemp
real :: A, B, C
integer :: i, j, k, n
real :: r

! saturation vapor pressure over water (C-C)

  if (do_CC) then

     do k=1,size(tabs,3)
        do j=1,size(tabs,2)
           do i=1,size(tabs,1)
              r = max(0.,min(2*num-1., (tabs(i,j,k) - T1)/dtemp))
              n = r
              r = r-n
              n = n-num
              esat(i,j,k) = (1.-r)*etable(n) + r*etable(n+1)
           enddo
        enddo
     enddo

     esat = exp(esat)
     return

  endif

! saturation vapor pressure over water (Tetens)

!  esat = atet*exp(btet*(tabs - 273.0)/(tabs - 36.0)))

  logtemp = log10(tabs)

  do k=1,size(tabs,3)
     do j=1,size(tabs,2)
        do i=1,size(tabs,1)

! saturation vapor pressure over water (Goff and Gratch)

           C  = tboil/tabs(i,j,k)
           E1 = D1*(C - 1.0)/C * ln10
           E2 = D2*(C - 1.0)   * ln10
           F1 = C5*(C - 1.0)
           F3 = C7*(exp(E1) - 1.0)
           F4 = C8*(exp(E2) - 1.0)
           F2 = C6*(logtb - logtemp(i,j,k))  ! c6*log10(c)
           B  = ln10*(F1 + F2 + F3 + F4 + F5)
           A  = 1.0

           esat(i,j,k) = A*exp(B)

           if (.not. present(wice)) cycle

! saturation vapor pressure over ice (Goff and Gratch)

           C  = tfreeze/tabs(i,j,k)
           F1 = C1*(C - 1.0)
           F3 = C3*(C - 1.0)/C
           F4 = C4
           F2 = C2*(logtf - logtemp(i,j,k))  ! c2*log10(c)
           B  = ln10*(F1 + F2 + F3 + F4)
           A  = 1.0

           esat(i,j,k) = esat(i,j,k)                                   & 
           + (A*exp(B) - esat(i,j,k))*wice(i,j,k)

        enddo
     enddo
  enddo

  return
end subroutine get_esat_3d

!#######################################################################

subroutine get_esat_2d ( tabs, esat, wice )

real, dimension(:,:), intent (in)  :: tabs
real, dimension(:,:), optional, intent (in)  :: wice
real, dimension(:,:), intent (out) :: esat

real, dimension(size(tabs,1),size(tabs,2)) :: logtemp
real :: A, B, C
integer :: j, k, n
real :: r

! saturation vapor pressure over water (C-C)

  if (do_CC) then
   
     do k=1,size(tabs,2)
        do j=1,size(tabs,1)
           r = max(0.,min(2*num-1., (tabs(j,k) - T1)/dtemp))
           n = r
           r = r-n
           n = n-num
           esat(j,k) = (1.-r)*etable(n) + r*etable(n+1)
        enddo
     enddo

     esat = exp(esat)
     return

  endif

! saturation vapor pressure over water (Tetens)

!  esat = atet*exp(btet*(tabs - 273.0)/(tabs - 36.0)))

  logtemp = log10(tabs)

  do k=1,size(tabs,2)
     do j=1,size(tabs,1)

! saturation vapor pressure over water (Goff and Gratch)

        C  = tboil/tabs(j,k)
        E1 = D1*(C - 1.0)/C * ln10
        E2 = D2*(C - 1.0)   * ln10
        F1 = C5*(C - 1.0)
        F3 = C7*(exp(E1) - 1.0)
        F4 = C8*(exp(E2) - 1.0)
        F2 = C6*(logtb - logtemp(j,k))  ! c6*log10(c)
        B  = ln10*(F1 + F2 + F3 + F4 + F5)
        A  = 1.0

        esat(j,k) = A*exp(B)

        if (.not. present(wice)) cycle

! saturation vapor pressure over ice (Goff and Gratch)

        C  = tfreeze/tabs(j,k)
        F1 = C1*(C - 1.0)
        F3 = C3*(C - 1.0)/C
        F4 = C4
        F2 = C2*(logtf - logtemp(j,k))  ! c2*log10(c)
        B  = ln10*(F1 + F2 + F3 + F4)
        A  = 1.0

        esat(j,k) = esat(j,k)                                          & 
        + (A*exp(B) - esat(j,k))*wice(j,k)

     enddo
  enddo

  return
end subroutine get_esat_2d

!#######################################################################

subroutine get_esat_0d ( tabs, esat )

real, intent (in)  :: tabs
real, intent (out) :: esat

real :: logtemp
real :: A, B, C
integer :: n
real :: r

! saturation vapor pressure over water (C-C)

  if (do_CC) then
   
     r = max(0.,min(2*num-1., (tabs - T1)/dtemp))
     n = r
     r = r-n
     n = n-num
     esat = (1.-r)*etable(n) + r*etable(n+1)

     esat = exp(esat)
     return

  endif

! saturation vapor pressure over water (Tetens)

  logtemp = log10(tabs)

! saturation vapor pressure over water (Goff and Gratch)

  C  = tboil/tabs
  E1 = D1*(C - 1.0)/C * ln10
  E2 = D2*(C - 1.0)   * ln10
  F1 = C5*(C - 1.0)
  F3 = C7*(exp(E1) - 1.0)
  F4 = C8*(exp(E2) - 1.0)
  F2 = C6*(logtb - logtemp)  ! c6*log10(c)
  B  = ln10*(F1 + F2 + F3 + F4 + F5)
  A  = 1.0

  esat = A*exp(B)

  return
end subroutine get_esat_0d

!#######################################################################

subroutine get_dlnesatdt ( tabs, dlnesatdt, wice )

real, dimension(:,:,:), intent (in)            :: tabs
real, dimension(:,:,:), optional, intent (in)  :: wice
real, dimension(:,:,:), intent (out)           :: dlnesatdt

real, dimension(size(tabs,1),size(tabs,2),size(tabs,3)) :: rtemp
real :: A, B, C
integer :: i, j, k

! log derivative of vapor pressure over water (C-C)

  if (do_CC) then
     dlnesatdt = Hlv/(Rvgas*tabs**2)
     return
  endif

! log derivative of vapor pressure over water (Tetens)

!  dlnesatdt = (273.0 - 36.0)*btet/(tabs - 36.0)**2

! log derivative of vapor pressure over water (Goff and Gratch)

  rtemp = 1.0/tabs

  do k=1,size(tabs,3)
     do j=1,size(tabs,2)
        do i=1,size(tabs,1)

           C  = tboil*rtemp(i,j,k)
           E1 = D1*(C - 1.0)/C * ln10
           E2 = D2*(C - 1.0)   * ln10
           F1 = C5
           F2 = C6/(ln10*C)
           F3 = C7*exp(E1)*ln10*D1/C**2
           F4 = C8*exp(E2)*ln10*D2
           B  = ln10*(F1 + F2 + F3 + F4)
           A  = -C*rtemp(i,j,k)

           dlnesatdt(i,j,k) = B*A

           if (.not. present(wice)) cycle

! log derivative of vapor pressure over ice (Goff and Gratch)

           C = tfreeze*rtemp(i,j,k)
           B = ln10*(C1 + (C2/ln10 + C3/C)/C)
           A = -C*rtemp(i,j,k)

           dlnesatdt(i,j,k) = dlnesatdt(i,j,k)                         &
                     + (B*A - dlnesatdt(i,j,k))*wice(i,j,k)

        enddo
     enddo
  enddo

  return
end subroutine get_dlnesatdt

!#######################################################################

subroutine moisture_init

real :: temp
integer :: n

integer :: io, ierr

namelist /moisture_nml/ do_kessler, do_ras, qvmin, do_CC

!-----------------------------------------------------------------------
! read namelist
!-----------------------------------------------------------------------

  read(input_nml_file, nml=moisture_nml, iostat=io)
  ierr = check_nml_error(io, 'moisture_nml')

!-----------------------------------------------------------------------
! write version number to logfile
!-----------------------------------------------------------------------

  call write_version_number (version, tag)

  logtf = log10(tfreeze)
  logtb = log10(tboil)

  do_init = .false.

  if (.not. do_CC ) return

  allocate (etable(-num:num))
  dtemp = Trange/num

  etable(0) = eref
  do n=1,num
     temp = T0 + (n-.5)*dtemp
     etable(n) = etable(n-1) + Hlv/(Rvgas*temp**2)*dtemp
     temp = T0 - (n-.5)*dtemp
     etable(-n) = etable(-n+1) - Hlv/(Rvgas*temp**2)*dtemp
  enddo

  T1 = T0 - Trange

  return
end subroutine moisture_init

!#######################################################################

end module zetac_moisture_mod
