module zetac_cgrid_interp_mod

!-----------------------------------------------------------------------
! module information
!-----------------------------------------------------------------------

use fms_mod,                   only : write_version_number
use mpp_mod,    only : mpp_pe

use zetac_horiz_grid_type_mod, only : horiz_grid_type

!-----------------------------------------------------------------------
! public interfaces
!-----------------------------------------------------------------------

implicit none
private
public cgrid_interp_uv, cgrid_interp_mass

interface cgrid_interp_uv
  module procedure cgrid_interp_uv_2d
  module procedure cgrid_interp_uv_3d
end interface

interface cgrid_interp_mass
  module procedure cgrid_interp_mass_2d
  module procedure cgrid_interp_mass_3d
end interface

character(len=*), parameter :: module='zetac_cgrid_interp_mod'
logical :: do_init=.true.

integer :: ibh, ieh, ibhw, iehe
integer :: jbh, jeh, jbhs, jehn

!-----------------------------------------------------------------------
! version control information
!-----------------------------------------------------------------------

character(len=128)  :: version =  '$Id: zetac_cgrid_interp.f90,v 1.1.2.8 2004/05/21 15:46:01 ck Exp $'
character(len=128)  :: tag     =  '$Name: latest $'

contains

!#######################################################################

subroutine cgrid_interp_uv_2d ( Mgrid, datam, datau, datav )

type (horiz_grid_type)                                        ::       &
                                                                Mgrid
     
real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed),             &
                                               intent (in)    ::       &
                                                                datam

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed),             &
                                               intent (out)   ::       &
                                                         datau, datav
  
!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

integer :: i, j

!-----------------------------------------------------------------------
! initialize module
!-----------------------------------------------------------------------

  if (do_init) then
     call cgrid_interp_init ( Mgrid )
  endif

  do j=jbh,jeh
     do i=ibh,iehe
        datau(i,j) = 0.5*(datam(i,j) + datam(i+1,j))
     enddo
  enddo
  if (iehe < ieh) then
     do j=jbh,jeh
        datau(ieh,j) = datam(ieh,j)
     enddo
  endif
  do j=jbh,jehn
     do i=ibh,ieh
        datav(i,j) = 0.5*(datam(i,j) + datam(i,j+1))
     enddo
  enddo
  if (jehn < jeh) then
     do i=ibh,ieh
        datav(i,jeh) = datam(i,jeh)
     enddo
  endif

  return
end subroutine cgrid_interp_uv_2d

!#######################################################################

subroutine cgrid_interp_uv_3d ( Mgrid, iz, datam, datau, datav )

type (horiz_grid_type)                                        ::       &
                                                                Mgrid
								
integer, intent(in)                                           ::   iz
     
real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed, iz),         &
                                               intent (in)    ::       &
                                                                datam

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed, iz),         &
                                               intent (out)   ::       &
                                                         datau, datav
  
!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

integer :: i, j

!-----------------------------------------------------------------------
! initialize module
!-----------------------------------------------------------------------

  if (do_init) then
    call cgrid_interp_init ( Mgrid )
  endif

  do j=jbh,jeh
     do i=ibh,iehe
        datau(i,j,:) = 0.5*(datam(i,j,:) + datam(i+1,j,:))
     enddo
  enddo
  if (iehe < ieh) then
     do j=jbh,jeh
        datau(ieh,j,:) = datam(ieh,j,:)
     enddo
  endif
  do j=jbh,jehn
     do i=ibh,ieh
        datav(i,j,:) = 0.5*(datam(i,j,:) + datam(i,j+1,:))
     enddo
  enddo
  if (jehn < jeh) then
     do i=ibh,ieh
        datav(i,jeh,:) = datam(i,jeh,:)
     enddo
  endif

  return
end subroutine cgrid_interp_uv_3d

!#######################################################################

subroutine cgrid_interp_mass_2d ( Mgrid, datau, datav, uatm, vatm )

type (horiz_grid_type)                                        ::       &
                                                                Mgrid

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed),             &
                                            intent (in)       ::       &
                                                         datau, datav

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed),             &
                                               intent (out)   ::       &
                                                           uatm, vatm

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

integer :: i, j
  
  if (do_init) then
    call cgrid_interp_init ( Mgrid )
  endif

  do j=jbh,jeh
     do i=ibhw,ieh
        uatm(i,j) = 0.5*(datau(i-1,j) + datau(i,j))
     enddo
  enddo

  do j=jbhs,jeh
     do i=ibh,ieh
        vatm(i,j) = 0.5*(datav(i,j-1) + datav(i,j))
     enddo
  enddo

!-----------------------------------------------------------------------
! extrapolate west and south in case boundary is not cyclic
!-----------------------------------------------------------------------

  if (ibh < ibhw) then
     do j=jbh,jeh
        uatm(ibh,j) = datau(ibh,j)
     enddo
  endif
  if (jbh < jbhs) then
     do i=ibh,ieh
        vatm(i,jbh) = datav(i,jbh)
     enddo
  endif

  return
end subroutine cgrid_interp_mass_2d

!#######################################################################

subroutine cgrid_interp_mass_3d ( Mgrid, iz, datau, datav, uatm, vatm )

type (horiz_grid_type)                                        ::       &
                                                                Mgrid
								
integer, intent(in)                                           ::   iz

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed, iz),         &
                                             intent (in)      ::       &
                                                         datau, datav

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed, iz),         &
                                             intent (out)     ::       &
                                                           uatm, vatm

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

integer :: i, j
  
  if (do_init) then
     call cgrid_interp_init ( Mgrid )
  endif

  do j=jbh,jeh
     do i=ibhw,ieh
        uatm(i,j,:) = 0.5*(datau(i-1,j,:) + datau(i,j,:))
     enddo
  enddo

  do j=jbhs,jeh
     do i=ibh,ieh
        vatm(i,j,:) = 0.5*(datav(i,j-1,:) + datav(i,j,:))
     enddo
  enddo

!-----------------------------------------------------------------------
! extrapolate west and south in case boundary is not cyclic
!-----------------------------------------------------------------------

  if ( ibhw > ibh ) then
     do j=jbh,jeh
        uatm(ibh,j,:) = datau(ibh,j,:)
     enddo
  endif
  if ( jbhs > jbh ) then
     do i=ibh,ieh
        vatm(i,jbh,:) = datav(i,jbh,:)
     enddo
  endif

  return
end subroutine cgrid_interp_mass_3d

!#######################################################################

subroutine cgrid_interp_init ( Mgrid )

type (horiz_grid_type)                                        ::       &
                                                                Mgrid

!-----------------------------------------------------------------------
! write version number to logfile
!-----------------------------------------------------------------------

  call write_version_number (version, tag)

  ibh  = Mgrid%ibh
  ieh  = Mgrid%ieh
  jbh  = Mgrid%jbh
  jeh  = Mgrid%jeh

  ibhw = max(ibh, Mgrid%ibg+1)
  jbhs = max(jbh, Mgrid%jbg+1)
  iehe = min(ieh, Mgrid%ieg-1)
  jehn = min(jeh, Mgrid%jeg-1)

  do_init = .false.
  
  return
end subroutine cgrid_interp_init

!#######################################################################

end module zetac_cgrid_interp_mod
