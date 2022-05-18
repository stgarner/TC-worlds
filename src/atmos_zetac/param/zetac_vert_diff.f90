module zetac_vert_diff_mod

!-----------------------------------------------------------------------
! module information
!-----------------------------------------------------------------------

use fms_mod,          only : write_version_number, mpp_pe, mpp_root_pe,&
                             check_nml_error, stdlog, error_mesg, FATAL
use mpp_mod,          only : input_nml_file

use constants_mod,    only : Grav

use field_manager_mod,           only : model_atmos

use zetac_horiz_grid_type_mod,   only : horiz_grid_type
use zetac_horiz_metric_type_mod, only : horiz_metric_type
use zetac_time_pointers_mod,     only : ntime
use zetac_tridiag_mod,           only : solve_tridiag
use zetac_update_halos_mod,      only : update_halos
use zetac_vert_metric_type_mod,  only : vert_metric_type

!-----------------------------------------------------------------------
! public interfaces
!-----------------------------------------------------------------------

implicit none
private
public vert_diff, vert_diff_init, vdiff_elements, vdiff_eval, vdifk_eval

character (len=*),  parameter :: module='zetac_vert_diff_mod'
character (len=16), parameter :: model='atmos_mod'

real, allocatable, dimension(:,:,:) :: diffu, diffl, du, dm, dl, dprat

! namelist:
real :: control=0.

integer :: ibc, iec
integer :: jbc, jec
integer :: kbd, ked

!-----------------------------------------------------------------------
! version control information
!-----------------------------------------------------------------------

character(len=128)  :: version =  '$Id: zetac_vert_diff.f90,v 1.1.2.9.2.4.2.5 2005/07/18 18:44:25 stg Exp $'
character(len=128)  :: tag     =  '$Name:  $'

contains

!#######################################################################

subroutine vert_diff (                                                 &
                                                                Mgrid, &
                                                             var, rhs, &
                                                   beta, gamma, ldiag )

type (horiz_grid_type),                        intent (in)    ::       &
                                                                Mgrid
       
real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),           intent (inout) ::      &
                                                                  var

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),           intent (in)    ::      &
                                                                  rhs

real, dimension(Mgrid%ibc:Mgrid%iec, Mgrid%jbc:Mgrid%jec,              &
                Mgrid%kbd+1:Mgrid%ked),         intent (inout) ::      &
                                                   beta, gamma, ldiag


  call solve_tridiag ( beta, gamma, ldiag,                             &
                       rhs = rhs(ibc:iec,jbc:jec,kbd+1:ked),           &
                       ans = var(ibc:iec,jbc:jec,kbd+1:ked) )

  return
end subroutine vert_diff

!#######################################################################

subroutine vdiff_elements (                                            &
                                                                Mgrid, &
                                                          diffw, nudg, &
                                                     dp, rfull, rhalf, &
                                                 cd, dflxdvar, gz, fz, &
                                                 rdzm, rdzw, dzm, dzw, &
                                                   beta, gamma, ldiag, &
                                             do_flux, present, future, &
                                                                 delt )

type (horiz_grid_type),                     intent (in)    ::          &
                                                                Mgrid
              
real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),                                  &
                                            intent (in)    ::          &
                                                          diffw, nudg

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked, ntime),                           &
                                            intent (in)    ::          &
                                                                   dp
real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),                                  &
                                            intent (in)    ::          &
                                                         rfull, rhalf

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed),             &
                                            intent (in)    ::          &
                                                         cd, dflxdvar

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),                                  &
                                            intent (in)    ::          &
                                                               gz, fz

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),                                  &
                                            intent (out)   ::          &
                                                           rdzm, rdzw

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),                                  &
                                            intent (out)   ::          &
                                                             dzm, dzw

real, dimension(Mgrid%ibc:Mgrid%iec, Mgrid%jbc:Mgrid%jec,              &
                Mgrid%kbd+1:Mgrid%ked),     intent (inout) ::          &
                                                   beta, gamma, ldiag

real,                                       intent (in)    ::          &
                                                                 delt

logical, intent (in) :: do_flux
integer, intent (in) :: present, future

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

integer :: i, j, k

!-----------------------------------------------------------------------
! inverse vertical grid size
!-----------------------------------------------------------------------

  do k=kbd+1,ked
     do j=jbc,jec
        do i=ibc,iec
           dzm(i,j,k) = (gz(i,j,k-1) - gz(i,j,k))/Grav
         rdzm(i,j,k) = 1./rfull(i,j,k)              ! divide rfull by dzm
        enddo
     enddo
  enddo

  do k=kbd+1,ked-1
     do j=jbc,jec
        do i=ibc,iec
           dzw(i,j,k) = (fz(i,j,k) - fz(i,j,k+1))/Grav
           rdzw(i,j,k) = rhalf(i,j,k)/dzw(i,j,k)**2 ! divide rhalf by dzw
        enddo
     enddo
  enddo
  do j=jbc,jec
     do i=ibc,iec
        dzw(i,j,ked) = (fz(i,j,ked) - gz(i,j,ked))/Grav
        rdzw(i,j,ked) = rhalf(i,j,ked)/dzw(i,j,ked) ! divide rhalf by dzw
     enddo
  enddo

  dzw(:,:,kbd) = dzw(:,:,kbd+1)
  rdzw(:,:,kbd) = rdzw(:,:,kbd+1)

!-----------------------------------------------------------------------
! tridiagonal elements
!-----------------------------------------------------------------------

  do k=kbd+1,ked
     do j=jbc,jec
        do i=ibc,iec
           diffu(i,j,k) = diffw(i,j,k)*rdzw(i,j,k)    *rdzm(i,j,k)
           diffl(i,j,k) = diffw(i,j,k-1)*rdzw(i,j,k-1)*rdzm(i,j,k)
           du(i,j,k) = -delt*diffu(i,j,k)
           dm(i,j,k) =  delt*(diffl(i,j,k) + diffu(i,j,k) + nudg(i,j,k))
           dl(i,j,k) = -delt*diffl(i,j,k)
        enddo
     enddo
  enddo

  do j=jbc,jec
     do i=ibc,iec
        dm(i,j,ked) = dm(i,j,ked)                                      &
              + delt*cd(i,j)*dflxdvar(i,j)*rdzm(i,j,ked)*rdzw(i,j,ked)
     enddo
  enddo

  if ( do_flux ) then
     dprat = dp(ibc:iec,jbc:jec,kbd+1:ked,present)/                    &
             dp(ibc:iec,jbc:jec,kbd+1:ked,future)
     du = du * dprat
     dm = dm * dprat
     dl = dl * dprat
  endif

  dm = 1.0 + dm

  ldiag = dl

  call solve_tridiag ( beta, gamma, dl, dm, du )

  return
end subroutine vdiff_elements

!#######################################################################

subroutine vdiff_eval (                                                &
                                                                Mgrid, &
                                                                diffw, &
                                                         cd, dflx, vs, &
                                                           rdzm, rdzw, &
                                                dvar, var, vbdy, tend, &
                                                         pres, future )

type (horiz_grid_type),                     intent (in)    ::          &
                                                                Mgrid
              
real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),                                  &
                                            intent (in)    ::          &
                                                                diffw

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed),             &
                                            intent (in)    ::          &
                                                         cd, dflx, vs

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),                                  &
                                            intent (in)    ::          &
                                                           rdzm, rdzw

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),                                  &
                                            intent (in)    ::          &
                                                                 dvar

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked,ntime), intent (in)    ::          &
                                                                  var

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed),             &
                                            intent (inout) ::          &
                                                                 vbdy

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),       intent (out)   ::          &
                                                                 tend

integer, intent (in) :: pres, future

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed)  ::          &
                                                                 vsfc
real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked)                       ::          &
                                                                 diff

integer :: i, j, k

!-----------------------------------------------------------------------
! diagnose diffusive tendency
!-----------------------------------------------------------------------

  do k=kbd+1,ked-1
     do j=jbc,jec
        do i=ibc,iec
           diff(i,j,k) = diffw(i,j,k)*rdzw(i,j,k)*                     &
                    (var(i,j,k,pres) - var(i,j,k+1,pres) + dvar(i,j,k))
        enddo
     enddo
  enddo

  diff(:,:,kbd) = 0.
  diff(:,:,ked) = 0.

  do k=kbd+1,ked
     do j=jbc,jec
        do i=ibc,iec
           tend(i,j,k) = (diff(i,j,k-1) - diff(i,j,k))*rdzm(i,j,k)
        enddo
     enddo
  enddo

  do j=jbc,jec
     do i=ibc,iec
        vbdy(i,j) = var(i,j,ked,future) - vbdy(i,j)
        tend(i,j,ked) = tend(i,j,ked) - cd(i,j)*                       &
                          vs(i,j)*vbdy(i,j)*rdzm(i,j,ked)*rdzw(i,j,ked)

     enddo
  enddo

end subroutine vdiff_eval

!#######################################################################

subroutine vdifk_eval (                                                &
                                                                Mgrid, &
                                                                diffw, &
                                                               cd, vs, &
                                                           rdzm, rdzw, &
                                                    var, tend1, tend2 )

type (horiz_grid_type),                     intent (in)    ::          &
                                                                Mgrid
              
real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),                                  &
                                            intent (in)    ::          &
                                                    diffw, rdzm, rdzw

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed),             &
                                            intent (in)    ::          &
                                                               cd, vs

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),       intent (in)    ::          &
                                                                  var

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),       intent (out)   ::          &
                                                         tend1, tend2

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked)                       ::          &
                                                                 diff

integer :: i, j, k

!-----------------------------------------------------------------------
! diagnose diffusive tendency
!-----------------------------------------------------------------------

  do k=kbd+1,ked-1
     do j=jbc,jec
        do i=ibc,iec
           diff(i,j,k) = 0.5*(var(i,j,k) + var(i,j,k+1))*             &
                diffw(i,j,k)*(var(i,j,k) - var(i,j,k+1))*rdzw(i,j,k)
        enddo
     enddo
  enddo

  diff(:,:,kbd) = 0.
  diff(:,:,ked) = 0.

  do k=kbd+1,ked
     do j=jbc,jec
        do i=ibc,iec
           tend1(i,j,k) = (diff(i,j,k-1) - diff(i,j,k))*rdzm(i,j,k)
        enddo
     enddo
  enddo

  do k=kbd+1,ked-1
     do j=jbc,jec
        do i=ibc,iec
           diff(i,j,k) = 0.5*(var(i,j,k) - var(i,j,k+1))*             &
                diffw(i,j,k)*(var(i,j,k) - var(i,j,k+1))*rdzw(i,j,k)
        enddo
     enddo
  enddo

  diff(:,:,kbd) = 0.
  diff(:,:,ked) = 0.

  do k=kbd+1,ked
     do j=jbc,jec
        do i=ibc,iec
           tend2(i,j,k) = -(diff(i,j,k-1) + diff(i,j,k))*rdzm(i,j,k)
        enddo
     enddo
  enddo

  do j=jbc,jec
     do i=ibc,iec
        tend2(i,j,ked) = tend2(i,j,ked) - cd(i,j)*                     &
                    vs(i,j)*var(i,j,ked)**2*rdzm(i,j,ked)*rdzw(i,j,ked)
     enddo
  enddo

end subroutine vdifk_eval

!#######################################################################

subroutine vert_diff_init ( Mgrid, Hmetric, Vmetric )

type (horiz_grid_type),                     intent (in)    ::          &
                                                                Mgrid
       
type (horiz_metric_type),                   intent (in)    ::          &
                                                              Hmetric
       
type (vert_metric_type),                    intent (in)    ::          &
                                                              Vmetric

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

integer :: io, ierr

namelist /vert_diff_nml/ &
  control

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

integer :: i, ibd, ied
integer :: j, jbd, jed

!-----------------------------------------------------------------------
! read namelist
!-----------------------------------------------------------------------

  read (input_nml_file, nml=vert_diff_nml, iostat=io)
  ierr = check_nml_error(io,'vert_diff_nml')

  ibc = Mgrid%ibc
  iec = Mgrid%iec
  jbc = Mgrid%jbc
  jec = Mgrid%jec

  ibd = Mgrid%ibd
  ied = Mgrid%ied
  jbd = Mgrid%jbd
  jed = Mgrid%jed
  kbd = Mgrid%kbd
  ked = Mgrid%ked

  allocate (diffu(ibc:iec,jbc:jec,kbd+1:ked))
  allocate (diffl(ibc:iec,jbc:jec,kbd+1:ked))
  allocate (   dl(ibc:iec,jbc:jec,kbd+1:ked))
  allocate (   dm(ibc:iec,jbc:jec,kbd+1:ked))
  allocate (   du(ibc:iec,jbc:jec,kbd+1:ked))
  allocate (dprat(ibc:iec,jbc:jec,kbd+1:ked))

!-----------------------------------------------------------------------
! write version number and namelist to logfile
!-----------------------------------------------------------------------

  call write_version_number (version, tag)

  if (mpp_pe() == mpp_root_pe())                                       &
                                    write (stdlog(), nml=vert_diff_nml)

  return
end subroutine vert_diff_init

!#######################################################################
   
subroutine error_handler ( message ) 
character(len=*), intent(in) :: message
   
   call error_mesg (module, message, FATAL)

end subroutine error_handler

!#######################################################################

end module zetac_vert_diff_mod
