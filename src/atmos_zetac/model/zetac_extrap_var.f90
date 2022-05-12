module zetac_extrap_var_mod

!-----------------------------------------------------------------------
! module information
!-----------------------------------------------------------------------

use mpp_mod,                    only : mpp_pe, mpp_root_pe
use fms_mod,                    only : write_version_number
use constants_mod,              only : grav, cp_air

use zetac_convert_var_mod,      only : get_tv
use zetac_horiz_grid_type_mod , only : horiz_grid_type
use zetac_vert_metric_type_mod, only : vert_metric_type

!-----------------------------------------------------------------------
! public interfaces
!-----------------------------------------------------------------------

implicit none
private
public extrap_ptq, extrap_uv, extrap_tracer, extrap_var_init

character(len=*), parameter :: module='zetac_extrap_var_mod'
logical :: do_lateral

integer :: ibc, iec, ibh, ieh, ibu
integer :: jbc, jec, jbh, jeh, jbv
integer :: kbd, ked
integer :: ix, iy

real :: buffac, fac, fac1
logical :: no_slip_w, no_slip_e, no_slip_s, no_slip_n

!-----------------------------------------------------------------------
! version control information
!-----------------------------------------------------------------------

character(len=128)  :: version =  '$Id: zetac_extrap_var.f90,v 1.1.2.6.2.6 2005/07/30 02:40:00 stg Exp $'
character(len=128)  :: tag     =  '$Name:  $'

contains

!#######################################################################

subroutine extrap_ptq (                                                &
                                                                Mgrid, &
                                                           ps, ta, qv, &
                                                  ppref, taref, qvref, &
                                                                 delt )

type (horiz_grid_type),                           intent (in)    ::    &
                                                                Mgrid

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed),             &
                                                  intent (inout) ::    &
                                                                   ps

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),             intent (inout) ::    &
                                                               ta, qv

real, dimension(Mgrid%jbd:Mgrid%jed),                                  &
                                                  intent (in)    ::    &
                                                                ppref
					   
real, dimension(Mgrid%jbd:Mgrid%jed, Mgrid%iz),                        &
                                                  intent (in)    ::    &
                                                         taref, qvref
					   
real,                                              intent (in)   ::    &
                                                                 delt

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

integer :: i, j, k

!-----------------------------------------------------------------------
! extrapolate to dummy mass level (use hydrostatic balance for pressure)
!-----------------------------------------------------------------------

  do j=jbc,jec
     ta(:,j,kbd) = taref(j,kbd)
     qv(:,j,kbd) = qvref(j,kbd)
  enddo

  if ( .not. do_lateral ) return

!-----------------------------------------------------------------------
! extrapolate to western and eastern halos
!-----------------------------------------------------------------------

  fac = min(1.0, delt*buffac)
  fac1 = 1.0 - fac

  if ( ibh < ibc ) then
     do i=ibh,ibc-1
        do j=jbc,jec
           ps(i,j) = ps(ibc,j)*fac1 + ppref(j)*fac
        enddo
     enddo
     do k=kbd,ked
        do i=ibh,ibc-1
	   do j=jbc,jec
              ta(i,j,k) = ta(ibc,j,k)*fac1 + taref(j,k)*fac
              qv(i,j,k) = qv(ibc,j,k)*fac1 + qvref(j,k)*fac
           enddo
        enddo
     enddo
  endif
  if ( ieh > iec ) then
     do i=iec+1,ieh
        do j=jbc,jec
           ps(i,j) = ps(iec,j)*fac1 + ppref(j)*fac
        enddo
     enddo
     do k=kbd,ked
        do i=iec+1,ieh
	   do j=jbc,jec
              ta(i,j,k) = ta(iec,j,k)*fac1 + taref(j,k)*fac
              qv(i,j,k) = qv(iec,j,k)*fac1 + qvref(j,k)*fac
           enddo
        enddo
     enddo
  endif

!-----------------------------------------------------------------------
! extrapolate to southern and northern halos
!-----------------------------------------------------------------------

  if ( jbh < jbc ) then
     do j=jbh,jbc-1
        do i=ibh,ieh
           ps(i,j) = ps(i,jbc)*fac1 + ppref(j)*fac
        enddo
     enddo
     do k=kbd,ked
        do j=jbh,jbc-1
           do i=ibh,ieh
              ta(i,j,k) = ta(i,jbc,k)*fac1 + taref(j,k)*fac
              qv(i,j,k) = qv(i,jbc,k)*fac1 + qvref(j,k)*fac
           enddo
	enddo
     enddo
  endif
  if ( jeh > jec ) then
     do j=jec+1,jeh
        do i=ibh,ieh
           ps(i,j) = ps(i,jec)*fac1 + ppref(j)*fac
        enddo
     enddo
     do k=kbd,ked
        do j=jec+1,jeh
           do i=ibh,ieh
              ta(i,j,k) = ta(i,jec,k)*fac1 + taref(j,k)*fac
              qv(i,j,k) = qv(i,jec,k)*fac1 + qvref(j,k)*fac
           enddo
        enddo
     enddo
  endif

  return
end subroutine extrap_ptq

!#######################################################################

subroutine extrap_uv ( Mgrid, uu, vv, uuref )

type (horiz_grid_type),                     intent (in)     ::         &
                                                                Mgrid

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),       intent (inout)  ::         &
                                                               uu, vv

real, dimension(Mgrid%jbd:Mgrid%jed, Mgrid%iz),                        &
                                            intent (in)     ::         &
                                                                uuref

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

integer :: i, j, k

!-----------------------------------------------------------------------
! extrapolate to dummy mass level
!-----------------------------------------------------------------------

  do j=jbc,jec
     uu(:,j,kbd) = uuref(j,kbd+1)
  enddo
  do j=jbv,jec
     vv(:,j,kbd) = 0.
  enddo

  if ( .not. do_lateral ) return

!-----------------------------------------------------------------------
! extrapolate to global halos
!-----------------------------------------------------------------------

  if ( ibh < ibc ) then
     do k=kbd,ked
        do i=ibh,ibu-1
           do j=jbc,jec
              uu(i,j,k) = uu(ibu,j,k)
           enddo
        enddo
        do i=ibh,ibc-1
           do j=jbv,jec
              vv(i,j,k) = vv(ibc,j,k)
           enddo
        enddo
     enddo
     if ( no_slip_w ) vv(ibh:ibc-1,jbv:jec,kbd:ked) = 0.0
  endif
  if ( ieh > iec ) then
     do k=kbd,ked
        do i=iec+1,ieh
           do j=jbc,jec
              uu(i,j,k) = uu(iec,j,k)
           enddo
        enddo
        do i=iec+1,ieh
           do j=jbv,jec
              vv(i,j,k) = vv(iec,j,k)
           enddo
        enddo
     enddo
     if ( no_slip_e ) vv(iec+1:ieh,jbv:jec,kbd:ked) = 0.0
  endif
  if ( jbh < jbc ) then
     do k=kbd,ked
        do j=jbh,jbc-1
           do i=ibh,ieh
              uu(i,j,k) = uu(i,jbc,k) + (uuref(j,k) - uuref(jbc,k))
           enddo
        enddo
        do j=jbh,jbv-1
           do i=ibh,ieh
              vv(i,j,k) = vv(i,jbv,k)
           enddo
        enddo
     enddo
     if ( no_slip_s ) uu(ibh:ieh,jbh:jbc-1,kbd:ked) = 0.0
  endif
  if ( jeh > jec ) then
     do k=kbd,ked
        do j=jec+1,jeh
           do i=ibh,ieh
              uu(i,j,k) = uu(i,jec,k) + (uuref(j,k) - uuref(jec,k))
              vv(i,j,k) = vv(i,jec,k)
           enddo
        enddo
     enddo
     if ( no_slip_n ) uu(ibh:ieh,jec+1:jeh,kbd:ked) = 0.0
  endif

  return
end subroutine extrap_uv

!#######################################################################

subroutine extrap_tracer ( Mgrid, qc, delt )

type (horiz_grid_type),                           intent (in)    ::    &
                                                                Mgrid

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),             intent (inout) ::    &
                                                                   qc

real, optional,                                    intent (in)   ::    &
                                                                 delt

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

integer :: i, j, k

!-----------------------------------------------------------------------
! extrapolate to dummy mass level
!-----------------------------------------------------------------------

  do j=jbc,jec
     qc(:,j,kbd) = 0.
  enddo

  if ( .not. do_lateral ) return

!-----------------------------------------------------------------------
! extrapolate to western and eastern halos
!-----------------------------------------------------------------------

  if (present(delt)) then
     fac = min(1.0, delt*buffac)
  else
     fac = 0.0
  endif
  fac1 = 1.0 - fac

  if ( ibh < ibc ) then
     do k=kbd,ked
        do i=ibh,ibc-1
	   do j=jbc,jec
              qc(i,j,k) = qc(ibc,j,k)*fac1
           enddo
        enddo
     enddo
  endif
  if ( ieh > iec ) then
     do k=kbd,ked
        do i=iec+1,ieh
	   do j=jbc,jec
              qc(i,j,k) = qc(iec,j,k)*fac1
           enddo
        enddo
     enddo
  endif

!-----------------------------------------------------------------------
! extrapolate to southern and northern halos
!-----------------------------------------------------------------------

  if ( jbh < jbc ) then
     do k=kbd,ked
        do j=jbh,jbc-1
           do i=ibh,ieh
              qc(i,j,k) = qc(i,jbc,k)*fac1
           enddo
	enddo
     enddo
  endif
  if ( jeh > jec ) then
     do k=kbd,ked
        do j=jec+1,jeh
           do i=ibh,ieh
              qc(i,j,k) = qc(i,jec,k)*fac1
           enddo
        enddo
     enddo
  endif

  return
end subroutine extrap_tracer

!#######################################################################

subroutine extrap_var_init (                                           &
                                             Mgrid, Vmetric, nudgbuff, &
                                       slip_w, slip_e, slip_s, slip_n, &
                                                               lhopen )

type (horiz_grid_type),                     intent (in)     ::         &
                                                                Mgrid

type (vert_metric_type),                    intent (in)     ::         &
                                                              Vmetric

real,                                       intent (in)     ::         &
                                                             nudgbuff

logical,                                    intent (in)     ::         &
                                       slip_w, slip_e, slip_s, slip_n

logical,                                    intent (in)     ::         &
                                                              lhopen

!-----------------------------------------------------------------------
! write version number to logfile
!-----------------------------------------------------------------------

  call write_version_number (version, tag)

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

  kbd = Mgrid%kbd
  ked = Mgrid%ked

  ix = Mgrid%ix
  iy = Mgrid%iy

  buffac = nudgbuff
  do_lateral = lhopen

  no_slip_w = .not. slip_w
  no_slip_e = .not. slip_e
  no_slip_s = .not. slip_s
  no_slip_n = .not. slip_n

  return
end subroutine extrap_var_init

!#######################################################################

end module zetac_extrap_var_mod
