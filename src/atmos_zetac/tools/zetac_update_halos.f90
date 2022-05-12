module zetac_update_halos_mod

!-----------------------------------------------------------------------
! module information
!-----------------------------------------------------------------------

use mpp_mod,                     only : mpp_send, mpp_recv,            &
                                        mpp_pe, mpp_root_pe, mpp_npes, &
                                        mpp_sync_self
use fms_mod,                     only : write_version_number
use mpp_domains_mod,             only : mpp_update_domains,            &
                                        mpp_get_layout, domain2D

use zetac_horiz_grid_type_mod,   only : horiz_grid_type

!-----------------------------------------------------------------------
! public interfaces
!-----------------------------------------------------------------------

implicit none

public update_halos

interface update_halos
  module procedure update_all
  module procedure update_3d
  module procedure update_2d
end interface

character(len=*), parameter :: module='zetac_update_halos_mod'
logical :: do_corners

integer :: ibc, iec, ibd, ied, ibh, ieh
integer :: jbc, jec, jbd, jed, jbh, jeh
integer :: nbuf
integer :: pe, pe_south, pe_north, pe_west, pe_east
integer :: npes_x, npes_y

type (domain2D), save :: Domain

!-----------------------------------------------------------------------
! version control information
!-----------------------------------------------------------------------

character(len=128)  :: version =  '$Id: zetac_update_halos.f90,v 1.1.2.1.2.7 2005/08/07 00:07:03 stg Exp $'
character(len=128)  :: tag     =  '$Name:  $'

contains

!#######################################################################

subroutine update_all (                                                &
                                                                Mgrid, &
                                               ps, uu, vv, ta, qv, qc, &
                                                               future )

type (horiz_grid_type), intent(in) :: Mgrid

real, dimension(:,:,:), intent (inout)   :: ps
real, dimension(:,:,:,:), intent (inout) :: uu, vv, ta, qv, qc

integer, intent(in) :: future

  call update_2d ( Mgrid, ps(:,:,future) )
  call update_3d ( Mgrid, uu(:,:,:,future) )
  call update_3d ( Mgrid, vv(:,:,:,future) )
  call update_3d ( Mgrid, ta(:,:,:,future) )
  call update_3d ( Mgrid, qv(:,:,:,future) )
  call update_3d ( Mgrid, qc(:,:,:,future) )

  return
end subroutine update_all

!#######################################################################

subroutine update_2d ( Mgrid, data )

type (horiz_grid_type), intent(in) :: Mgrid

real, dimension(:,:), intent (inout) :: data

!-----------------------------------------------------------------------
! update domain boundaries
!-----------------------------------------------------------------------

  call mpp_update_domains ( data, Domain )

!-----------------------------------------------------------------------
! fill in data domain corners at open boundaries
!-----------------------------------------------------------------------

  if ( do_corners ) call update_corners_2d ( Mgrid, data )

  return
end subroutine update_2d

!#######################################################################

subroutine update_3d ( Mgrid, data )

type (horiz_grid_type), intent(in) :: Mgrid

real, dimension(:,:,:), intent (inout) :: data

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

integer :: nz

!-----------------------------------------------------------------------
! update domain boundaries
!-----------------------------------------------------------------------

  call mpp_update_domains ( data, Domain )

!-----------------------------------------------------------------------
! fill in data domain corners at open boundaries
!-----------------------------------------------------------------------

  if ( do_corners ) then
     nz = size(data,3)
     call update_corners_3d ( Mgrid, data, nz )
  endif

  return
end subroutine update_3d 

!#######################################################################

subroutine update_corners_2d ( Mgrid, data )

type (horiz_grid_type), intent(in) :: Mgrid

real, dimension(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed),              &
                                                intent (inout) :: data

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

real, dimension(nbuf, nbuf) :: sw_send, nw_recv
real, dimension(nbuf, nbuf) :: nw_send, sw_recv
real, dimension(nbuf, nbuf) :: se_send, ne_recv
real, dimension(nbuf, nbuf) :: ne_send, se_recv
real, dimension(nbuf, nbuf) :: ws_send, es_recv
real, dimension(nbuf, nbuf) :: es_send, ws_recv
real, dimension(nbuf, nbuf) :: wn_send, en_recv
real, dimension(nbuf, nbuf) :: en_send, wn_recv

integer :: length

!-----------------------------------------------------------------------
! update corners at open boundaries
!-----------------------------------------------------------------------

  length = nbuf*nbuf

  if ( ibh < ibc ) then

     if ( npes_y == 1 .and. .not. Mgrid%lyopen ) then
        data(ibd:ibc-1,jec+1:jed) = data(ibd:ibc-1,jbc:jbc+nbuf-1)
        data(ibd:ibc-1,jbd:jbc-1) = data(ibd:ibc-1,jec-nbuf+1:jec)
     else
        if ( jbc == jbh ) then
	   sw_send = data(ibd:ibc-1,jbc:jbc+nbuf-1)
           call mpp_send ( sw_send, length, pe_south )
        endif
        if ( jec == jeh ) then
           call mpp_recv ( nw_recv, length, pe_north )
           data(ibd:iec-1,jec+1:jed) = nw_recv
        endif
        if ( jec == jeh ) then
           nw_send = data(ibd:iec-1,jec-nbuf+1:jec)
           call mpp_send ( nw_send, length, pe_north )
        endif
        if ( jbc == jbh ) then
           call mpp_recv ( sw_recv, length, pe_south )
           data(ibd:iec-1,jbd:jbc-1) = sw_recv
        endif
     endif

  endif

  if ( ieh > iec ) then

     if ( npes_y == 1 .and. .not. Mgrid%lyopen ) then
        data(iec+1:ied,jec+1:jed) = data(iec+1:ied,jbc:jbc+nbuf-1)
        data(iec+1:ied,jbd:jbc-1) = data(iec+1:ied,jec-nbuf+1:jec)
     else
        if ( jbc == jbh ) then
           se_send = data(iec+1:ied,jbc:jbc+nbuf-1)
           call mpp_send ( se_send, length, pe_south )   
        endif
        if ( jec == jeh ) then
           call mpp_recv ( ne_recv, length, pe_north )
           data(iec+1:ied,jec+1:jed) = ne_recv
        endif
        if ( jec == jeh ) then
           ne_send = data(iec+1:ied,jec-nbuf+1:jec)
           call mpp_send ( ne_send, length, pe_north )
        endif
        if ( jbc == jbh ) then
           call mpp_recv ( se_recv, length, pe_south )
           data(iec+1:ied,jbd:jbc-1) = se_recv
        endif
     endif

  endif

  if ( jbh < jbc ) then

     if ( npes_x == 1 .and. .not. Mgrid%lxopen ) then
        data(iec+1:ied,jbd:jbc-1) = data(ibc:ibc+nbuf-1,jbd:jbc-1)
        data(ibd:ibc-1,jbd:jbc-1) = data(iec-nbuf+1:iec,jbd:jbc-1)
     else
        if ( ibc == ibh ) then
           ws_send = data(ibc:ibc+nbuf-1,jbd:jbc-1)
           call mpp_send ( ws_send, length, pe_west )
        endif
        if ( iec == ieh ) then
           call mpp_recv ( es_recv, length, pe_east )
           data(iec+1:ied,jbd:jbc-1) = es_recv
        endif
        if ( iec == ieh ) then
           es_send = data(iec-nbuf+1:iec,jbd:jbc-1)
           call mpp_send ( es_send, length, pe_east )
        endif
        if ( ibc == ibh ) then
           call mpp_recv ( ws_recv, length, pe_west )
           data(ibd:ibc-1,jbd:jbc-1) = ws_recv
        endif
     endif

  endif

  if ( jeh > jec ) then

     if ( npes_x == 1 .and. .not. Mgrid%lxopen ) then
        data(iec+1:ied,jec+1:jed) = data(ibc:ibc+nbuf-1,jec+1:jed)
        data(ibd:ibc-1,jec+1:jed) = data(iec-nbuf+1:iec,jec+1:jed)
     else
        if ( ibc == ibh ) then
           wn_send = data(ibc:ibc+nbuf-1,jec+1:jed)
           call mpp_send ( wn_send, length, pe_west )
        endif
        if ( iec == ieh ) then
           call mpp_recv ( en_recv, length, pe_east )
           data(iec+1:ied,jec+1:jed) = en_recv
        endif
        if ( iec == ieh ) then
           en_send = data(iec-nbuf+1:iec,jec+1:jed)
           call mpp_send ( en_send, length, pe_east )
        endif
        if ( ibc == ibh ) then
           call mpp_recv ( wn_recv, length, pe_west )
           data(ibd:ibc-1,jec+1:jed) = wn_recv
        endif
     endif

  endif

  call mpp_sync_self

  return
end subroutine update_corners_2d

!#######################################################################

subroutine update_corners_3d ( Mgrid, data, nz )

type (horiz_grid_type), intent(in) :: Mgrid

real, dimension(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed,nz),           &
                                                intent (inout) :: data

integer, intent( in) :: nz

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

real, dimension(nbuf, nbuf, nz) :: sw_send, nw_recv
real, dimension(nbuf, nbuf, nz) :: nw_send, sw_recv
real, dimension(nbuf, nbuf, nz) :: se_send, ne_recv
real, dimension(nbuf, nbuf, nz) :: ne_send, se_recv
real, dimension(nbuf, nbuf, nz) :: ws_send, es_recv
real, dimension(nbuf, nbuf, nz) :: es_send, ws_recv
real, dimension(nbuf, nbuf, nz) :: wn_send, en_recv
real, dimension(nbuf, nbuf, nz) :: en_send, wn_recv

integer :: length

!-----------------------------------------------------------------------
! update corners at open boundaries
!-----------------------------------------------------------------------

  length = nbuf*nbuf*nz

  if ( ibh < iec ) then

     if ( npes_y == 1 .and. .not. Mgrid%lyopen ) then
        data(ibd:ibc-1,jec+1:jed,:) = data(ibd:ibc-1,jbc:jbc+nbuf-1,:)
        data(ibd:ibc-1,jbd:jbc-1,:) = data(ibd:ibc-1,jec-nbuf+1:jec,:)
     else
        if ( jbc == jbh ) then
           sw_send = data(ibd:iec-1,jbc:jbc+nbuf-1,:)
           call mpp_send ( sw_send, length, pe_south )
        endif
        if ( jec == jeh ) then
           call mpp_recv ( nw_recv, length, pe_north )
           data(ibd:iec-1,jec+1:jed,:) = nw_recv
        endif
        if ( jec == jeh ) then
           nw_send = data(ibd:iec-1,jec-nbuf+1:jec,:)
           call mpp_send ( nw_send, length, pe_north )
        endif
        if ( jbc == jbh ) then
           call mpp_recv ( sw_recv, length, pe_south )
           data(ibd:iec-1,jbd:jbc-1,:) = sw_recv
        endif
     endif

  endif

  if ( ieh > iec ) then

     if ( npes_y == 1 .and. .not. Mgrid%lyopen ) then
        data(iec+1:ied,jec+1:jed,:) = data(iec+1:ied,jbc:jbc+nbuf-1,:)
        data(iec+1:ied,jbd:jbc-1,:) = data(iec+1:ied,jec-nbuf+1:jec,:)
     else
        if ( jbc == jbh ) then
           se_send = data(iec+1:ied,jbc:jbc+nbuf-1,:)
           call mpp_send ( se_send, length, pe_south )
        endif
        if ( jec == jeh ) then
           call mpp_recv ( ne_recv, length, pe_north )
           data(iec+1:ied,jec+1:jed,:) = ne_recv
        endif
        if ( jec == jeh ) then
           ne_send = data(iec+1:ied,jec-nbuf+1:jec,:)
           call mpp_send ( ne_send, length, pe_north )
        endif
        if ( jbc == jbh ) then
           call mpp_recv ( se_recv, length, pe_south )
           data(iec+1:ied,jbd:jbc-1,:) = se_recv
        endif
     endif

  endif

  if ( jbh < jbc ) then

     if ( npes_x == 1 .and. .not. Mgrid%lxopen ) then
        data(iec+1:ied,jbd:jbc-1,:) = data(ibc:ibc+nbuf-1,jbd:jbc-1,:)
        data(ibd:ibc-1,jbd:jbc-1,:) = data(iec-nbuf+1:iec,jbd:jbc-1,:)
     else
        if ( ibc == ibh ) then
           ws_send = data(ibc:ibc+nbuf-1,jbd:jbc-1,:)
           call mpp_send ( ws_send, length, pe_west )
        endif
        if ( iec == ieh ) then
           call mpp_recv ( es_recv, length, pe_east )
           data(iec+1:ied,jbd:jbc-1,:) = es_recv
        endif
        if ( iec == ieh ) then
           es_send = data(iec-nbuf+1:iec,jbd:jbc-1,:)
           call mpp_send ( es_send, length, pe_east )
        endif
        if ( ibc == ibh ) then
           call mpp_recv ( ws_recv, length, pe_west )
           data(ibd:ibc-1,jbd:jbc-1,:) = ws_recv
        endif
     endif

  endif

  if ( jeh > jec ) then

     if ( npes_x == 1 .and. .not. Mgrid%lxopen ) then
        data(iec+1:ied,jec+1:jed,:) = data(ibc:ibc+nbuf-1,jec+1:jed,:)
        data(ibd:ibc-1,jec+1:jed,:) = data(iec-nbuf+1:iec,jec+1:jed,:)
     else
        if ( ibc == ibh ) then
           wn_send = data(ibc:ibc+nbuf-1,jec+1:jed,:)
           call mpp_send ( wn_send, length, pe_west )
        endif
        if ( iec == ieh ) then
           call mpp_recv ( en_recv, length, pe_east )
           data(iec+1:ied,jec+1:jed,:) = en_recv
	endif
        if ( iec == ieh ) then
           en_send = data(iec-nbuf+1:iec,jec+1:jed,:)
           call mpp_send ( en_send, length, pe_east )
        endif
        if ( ibc == ibh ) then
           call mpp_recv ( wn_recv, length, pe_west )
           data(ibd:ibc-1,jec+1:jed,:) = wn_recv
        endif
     endif

  endif

  call mpp_sync_self

  return
end subroutine update_corners_3d

!#######################################################################

subroutine update_halos_init ( Mgrid )

type (horiz_grid_type), intent(in) :: Mgrid

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

integer :: layout(2), global(4), npes

!-----------------------------------------------------------------------
! write version number to logfile
!-----------------------------------------------------------------------

  call write_version_number (version, tag)

  call mpp_get_layout ( Mgrid%Domain, layout )

  npes_x = layout(1)
  npes_y = layout(2)
  npes = npes_x*npes_y

  do_corners = ( npes > 1 )
  if ( .not. do_corners ) return

  nbuf = Mgrid%nbuf
  Domain = Mgrid%Domain

  do_corners = (Mgrid%lxopen .or. Mgrid%lyopen)

  if ( .not. do_corners ) return

  ibc = Mgrid%ibc
  iec = Mgrid%iec
  ibd = Mgrid%ibd
  ied = Mgrid%ied
  ibh = Mgrid%ibh
  ieh = Mgrid%ieh

  jbc = Mgrid%jbc
  jec = Mgrid%jec
  jbd = Mgrid%jbd
  jed = Mgrid%jed
  jbh = Mgrid%jbh
  jeh = Mgrid%jeh

  pe = mpp_pe()

  if ( mod(pe,npes_x) > 0 ) then
     pe_west = pe - 1
  else
     pe_west = pe - 1 + npes_x
  endif
  if ( mod(pe+1,npes_x) > 0 ) then
     pe_east = pe + 1
  else
     pe_east = pe + 1 - npes_x
  endif

  if ( pe >= npes_x ) then
     pe_south = pe - npes_x
  else
     pe_south = pe - npes_x + npes
  endif
  if ( pe + npes_x < npes ) then
     pe_north = pe + npes_x
  else
     pe_north = pe + npes_x - npes
  endif

  return
end subroutine update_halos_init

!#######################################################################

end module zetac_update_halos_mod
