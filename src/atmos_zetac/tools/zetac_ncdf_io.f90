module zetac_ncdf_io_mod

!-----------------------------------------------------------------------
! module information
!-----------------------------------------------------------------------

use mpp_mod,                   only : mpp_pe, mpp_npes,                &
                                      mpp_sum, mpp_root_pe
use mpp_io_mod,                only : fieldtype, mpp_read, mpp_write
use fms_mod,                   only : write_version_number

use mpp_domains_mod,           only : mpp_get_compute_domain,          &
                                      mpp_get_data_domain,             &
                                      mpp_update_domains, domain2D,    &
                                      mpp_global_field,                &
                                      xupdate, yupdate

use diag_manager_mod,          only : send_data, need_data
use time_manager_mod,          only : time_type, operator(+)

use zetac_axes_mod,            only : gridu, gridv, gridm
use zetac_horiz_grid_type_mod, only : horiz_grid_type

!-----------------------------------------------------------------------
! public interfaces
!-----------------------------------------------------------------------

implicit none
private
public ncdf_read, ncdf_write, ncdf_fms_write, ncdf_fms_write_integral, &
       ncdf_io_init

interface ncdf_fms_write
  module procedure ncdf_fms_write_2d
  module procedure ncdf_fms_write_3d
end interface

character(len=*), parameter :: module='zetac_ncdf_io_mod'

integer :: ibu, ieu, ibu1
integer :: jbv, jev, jbv1
integer :: ib1, ie1, ib2, ie2
integer :: jb1, je1, jb2, je2
integer :: ibc, iec, ibh, ieh, ibd, ied
integer :: jbc, jec, jbh, jeh, jbd, jed
integer :: kbd, ked

logical :: ucover, vcover

integer, allocatable, dimension(:) :: test

real, allocatable, dimension(:,:)   :: data2, data2u, data2v, data2m
real, allocatable, dimension(:,:,:) :: data3, data3u, data3v, data3m

type (time_type), save :: Time_step

!-----------------------------------------------------------------------
! version control information
!-----------------------------------------------------------------------

character(len=128)  :: version =  '$Id: zetac_ncdf_io.f90,v 1.1.2.9.2.4 2005/07/05 02:11:03 stg Exp $'
character(len=128)  :: tag     =  '$Name:  $'

contains

!#######################################################################

subroutine ncdf_read ( Mgrid, unit, Field, kb, ke, data_out, timelevel )

type (horiz_grid_type),                           intent (in)     ::   &
                                                                Mgrid

integer,                                          intent (in)     ::   &
                                                                 unit

type (fieldtype),                                 intent (in)     ::   &
                                                                Field

integer,                                           intent(in)     ::   &
                                                               kb, ke

real, dimension(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed,kb:ke),        &
                                                  intent (out)    ::   &
                                                             data_out

integer, optional,                                intent(in)      ::   &
                                                            timelevel

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

type (domain2D) :: Domain        

real, dimension(ib2:ie2,jb2:je2,kb:ke)                            ::   &
                                                                 data

integer :: i, j, k

  call mpp_read ( unit, Field, Mgrid%Domain_plusxy, data, timelevel )

  Domain = Mgrid%Domain_plusxy
  call mpp_update_domains ( data, Domain )

  do k=kb,ke
     do j=jbh,jeh
        do i=ibh,ieh
           data_out(i,j,k) = data(i,j,k)
        enddo
     enddo
  enddo

  return
end subroutine ncdf_read

!#######################################################################

subroutine ncdf_write ( Mgrid, unit, Field, kb, ke, data_in, time )

type (horiz_grid_type),                           intent (inout)  ::   &
                                                                Mgrid

integer,                                          intent (in)     ::   &
                                                                 unit

type (fieldtype),                                 intent (in)     ::   &
                                                                Field

integer,                                          intent (in)     ::   &
                                                               kb, ke

real, dimension(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed,kb:ke),        &
                                                  intent (in)     ::   &
                                                              data_in

real, optional,                                   intent (in)     ::   &
                                                                 time

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

real, dimension(ib1:ie1,jb1:je1,kb:ke)                            ::   &
                                                                 data

integer :: i, j, k
 
  do k=kb,ke
     do j=jb1,je1
        do i=ib1,ie1
           data(i,j,k) = data_in(i,j,k)
        enddo
     enddo
  enddo

  if ( present(time) ) then
     call mpp_write ( unit, Field, Mgrid%Domain_plusxy, data, time )
  else
     call mpp_write ( unit, Field, Mgrid%Domain_plusxy, data )
  endif

  return
end subroutine ncdf_write

!#######################################################################

subroutine ncdf_fms_write_2d ( Mgrid, Time, id_field, grid, data_in )

type (horiz_grid_type),                           intent (in)     ::   &
                                                                Mgrid

type (time_type),                                 intent (in)     ::   &
                                                                 Time

real, dimension(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed),              &
                                                  intent (in)     ::   &
                                                              data_in

integer,                                          intent (in)     ::   &
                                                       id_field, grid

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

logical :: used
integer :: i, j

  if ( id_field <= 0 ) return
  if ( .not. need_data(id_field, Time + Time_step) ) return

  if ( grid == gridu ) then

     if (ucover) then
     
        call mpp_global_field ( Mgrid%Domain, data_in, data2,          &
                                                        flags=xupdate )

        do j=jbc,jec
           do i=ibu1,ieu
              data2u(i,j) = data2(i,j)
           enddo
        enddo
        if ( ibu == 0 ) then
           do j=jbc,jec
              data2u(ibu,j) = data_in(ibu,j)
           enddo
        endif

     else

        do j=jbc,jec
           do i=ibu,ieu
              data2u(i,j) = data_in(i,j)
           enddo
        enddo

     endif
        
     used = send_data ( id_field, data2u, Time )

  else if ( grid == gridv ) then

     if (vcover) then

        call mpp_global_field ( Mgrid%Domain, data_in, data2,          &
                                                        flags=yupdate )

        do j=jbv1,jev
           do i=ibc,iec
              data2v(i,j) = data2(i,j)
           enddo
        enddo
        if ( jbv == 0 ) then
           do i=ibc,iec
              data2v(i,jbv) = data_in(i,jbv)
           enddo
        endif

     else
     
        do j=jbv,jev
           do i=ibc,iec
              data2v(i,j) = data_in(i,j)
           enddo
        enddo

     endif

     used = send_data ( id_field, data2v, Time )

  else

     do j=jbc,jec
        do i=ibc,iec
           data2m(i,j) = data_in(i,j)
        enddo
     enddo

     used = send_data ( id_field, data2m, Time )

  endif
  
  return
end subroutine ncdf_fms_write_2d

!#######################################################################

subroutine ncdf_fms_write_3d ( Mgrid, Time, id_field, grid, data_in )

type (horiz_grid_type),                           intent (in)     ::   &
                                                                Mgrid

type (time_type),                                 intent (in)     ::   &
                                                                 Time

real, dimension(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed,               &
                Mgrid%kbd:Mgrid%ked),                                  &
                                                  intent (in)     ::   &
                                                              data_in

integer,                                          intent (in)     ::   &
                                                       id_field, grid

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

logical :: used
integer :: i, j, k

  if ( id_field <= 0 ) return
  if ( .not. need_data(id_field, Time + Time_step) ) return

  if ( grid == gridu ) then

     if (ucover) then
     
        call mpp_global_field ( Mgrid%Domain, data_in, data3,          &
                                                        flags=xupdate )

        do k=kbd,ked
           do j=jbc,jec
              do i=ibu1,ieu
                 data3u(i,j,k) = data3(i,j,k)
              enddo
           enddo
        enddo
        if ( ibu == 0 ) then
           do k=kbd,ked
              do j=jbc,jec
                 data3u(ibu,j,k) = data_in(ibu,j,k)
              enddo
           enddo
        endif
     
     else
     
        do k=kbd,ked
           do j=jbc,jec
              do i=ibu,ieu
                 data3u(i,j,k) = data_in(i,j,k)
              enddo
           enddo
        enddo

     endif     

     used = send_data ( id_field, data3u, Time )

  else if ( grid == gridv ) then

     if (vcover) then

        call mpp_global_field ( Mgrid%Domain, data_in, data3,          &
                                                        flags=yupdate )

        do k=kbd,ked
           do j=jbv1,jev
              do i=ibc,iec
                 data3v(i,j,k) = data3(i,j,k)
	      enddo
           enddo
        enddo
        if ( jbv == 0 ) then
           do k=kbd,ked
              do i=ibc,iec
                 data3v(i,jbv,k) = data_in(i,jbv,k)
              enddo
           enddo
        endif

     else
     
        do k=kbd,ked
           do j=jbv,jev
              do i=ibc,iec
                 data3v(i,j,k) = data_in(i,j,k)
	      enddo
           enddo
        enddo

     endif

     used = send_data ( id_field, data3v, Time )

  else 

     do k=kbd,ked
        do j=jbc,jec
           do i=ibc,iec
              data3m(i,j,k) = data_in(i,j,k)
           enddo
        enddo
     enddo

     used = send_data ( id_field, data3m, Time )

  endif

  return
end subroutine ncdf_fms_write_3d

!#######################################################################

subroutine ncdf_fms_write_integral (                                   &
                                                          Mgrid, Time, &
                                                             id_field, &
                                            data_in, dzee, dens, gfac )

type (horiz_grid_type),                           intent (in)     ::   &
                                                                Mgrid

type (time_type),                                 intent (in)     ::   &
                                                                 Time

integer,                                          intent (in)     ::   &
                                                             id_field

real, dimension(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed,               &
                Mgrid%kbd:Mgrid%ked),                                  &
                                                  intent (in)     ::   &
                                                        data_in, dzee

real, dimension(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed,               &
                Mgrid%kbd:Mgrid%ked), optional,                        &
                                                  intent (in)     ::   &
                                                                 dens

real, optional,                                   intent (in)     ::   &
                                                                 gfac

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

real, dimension(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed)          ::   &
                                                                 data
integer :: i, j

  if ( id_field <= 0 ) return
  if ( .not. need_data(id_field, Time + Time_step) ) return

  if (present(dens)) then

     do j=jbc,jec
        do i=ibc,iec
           data(i,j) = sum ( data_in(i,j,kbd+1:ked)*                   &
                                dens(i,j,kbd+1:ked)*                   &
                                dzee(i,j,kbd+1:ked) )
        enddo
     enddo

  else

     do j=jbc,jec
        do i=ibc,iec
           data(i,j) = sum ( data_in(i,j,kbd+1:ked)*                   &
                                dzee(i,j,kbd+1:ked) )
        enddo
     enddo

  endif

  if ( present(gfac) ) data = data/gfac
  
  call ncdf_fms_write_2d ( Mgrid, Time, id_field, gridm, data )

  return
end subroutine ncdf_fms_write_integral

!#######################################################################

subroutine ncdf_io_init ( Mgrid, Time_step_in )

type (horiz_grid_type),                           intent (in)     ::   &
                                                                Mgrid

type (time_type),                                 intent (in)     ::   &
                                                         Time_step_in

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

integer :: pe, npe
integer :: xsize, xmax_size, ysize, ymax_size
integer :: xbegin, xend, ybegin, yend

  ibc = Mgrid%ibc
  iec = Mgrid%iec
  ibh = Mgrid%ibh
  ieh = Mgrid%ieh

  jbc = Mgrid%jbc
  jec = Mgrid%jec
  jbh = Mgrid%jbh
  jeh = Mgrid%jeh
  
  ibd = Mgrid%ibd
  ied = Mgrid%ied
  jbd = Mgrid%jbd
  jed = Mgrid%jed

  kbd = Mgrid%kbd
  ked = Mgrid%ked

!-----------------------------------------------------------------------
! dcomposition for restart fields
!-----------------------------------------------------------------------

  call mpp_get_compute_domain (                                        &
                                                  Mgrid%Domain_plusxy, &
                                           xbegin, xend, ybegin, yend, &
                                   xsize, xmax_size, ysize, ymax_size )

  ib1 = xbegin
  ie1 = xend
  jb1 = ybegin
  je1 = yend

  call mpp_get_data_domain (                                           &
                                                  Mgrid%Domain_plusxy, &
                                           xbegin, xend, ybegin, yend, &
                                   xsize, xmax_size, ysize, ymax_size )

  ib2 = xbegin
  ie2 = xend
  jb2 = ybegin
  je2 = yend

!-----------------------------------------------------------------------
! decompositions for history fields
!-----------------------------------------------------------------------

  call mpp_get_compute_domain (                                        &
                                                       Mgrid%Domain_u, &
                                           xbegin, xend, ybegin, yend, &
                                   xsize, xmax_size, ysize, ymax_size )

  ibu = xbegin
  ieu = xend
  ibu1 = max(1,ibu)

  call mpp_get_compute_domain (                                        &
                                                       Mgrid%Domain_v, &
                                           xbegin, xend, ybegin, yend, &
                                   xsize, xmax_size, ysize, ymax_size )

  jbv = ybegin
  jev = yend
  jbv1 = max(1,jbv)

!-----------------------------------------------------------------------
! allocate work arrays
!-----------------------------------------------------------------------

  allocate (data2(Mgrid%ix,Mgrid%iy))
  allocate (data3(Mgrid%ix,Mgrid%iy,Mgrid%iz))

  allocate (data3u(ibu:ieu,jbc:jec,kbd:ked))
  allocate (data3v(ibc:iec,jbv:jev,kbd:ked))
  allocate (data3m(ibc:iec,jbc:jec,kbd:ked))
  allocate (data2u(ibu:ieu,jbc:jec))
  allocate (data2v(ibc:iec,jbv:jev))
  allocate (data2m(ibc:iec,jbc:jec))

  pe = mpp_pe()
  npe = mpp_npes()

  allocate (test(0:npe))

  test = 0
  if ( ibu >= ibd .and. ieu <= ied ) test(pe) = 1
  call mpp_sum (test, length=npe)
  ucover = any(test == 0)

  test = 0
  if ( jbv >= jbd .and. jev <= jed ) test(pe) = 1
  call mpp_sum (test, length=npe)
  vcover = any(test == 0)
  
  deallocate (test)
  
!-----------------------------------------------------------------------
! write version number to logfile
!-----------------------------------------------------------------------

  call write_version_number (version, tag)

  Time_step = Time_step_in
 
  return
end subroutine ncdf_io_init

!#######################################################################

end module zetac_ncdf_io_mod
