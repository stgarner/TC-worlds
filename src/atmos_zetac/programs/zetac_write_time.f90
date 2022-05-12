program zetac_write_stuff

use mpp_mod, only :   mpp_init, mpp_exit, mpp_pe, mpp_root_pe

use fms_mod, only :   fms_init, file_exist, open_namelist_file,        &
                      write_version_number,                            &
                      check_nml_error, stdlog, stdout,                 &
                      open_namelist_file, open_file, close_file

use time_manager_mod, only : time_type, set_time, set_date

use zetac_namelist_mod, only : namelist_read
use zetac_times_mod,    only : get_times, write_timestamps

implicit none

type (time_type) :: Time_init, Time_end

integer :: ix, iy, iz, nsteps, nbuf, nseg
logical :: lxopen, lyopen, lcartesian, lcylindrical
real    :: rlonmin, rlonmax, rlatmin, rlatmax
real    :: pbot, ptop
real, dimension(20) :: zeta
real, dimension(:), allocatable :: zetaw

integer :: ierr, unit, io
integer :: num_segs=1

namelist /segments_nml/ num_segs

!-----------------------------------------------------------------------
! version control information
!-----------------------------------------------------------------------

character(len=128)  :: version =  '$Id: zetac_write_stuff.f90,v 1.1.2.1.2.3 2004/08/25 19:50:23 ck Exp $'
character(len=128)  :: tag     =  '$Name:  $'

!-----------------------------------------------------------------------
! initialize fms and read segments namelist
!-----------------------------------------------------------------------

  call fms_init ()

  if (file_exist('input.nml')) then

    unit =  open_namelist_file ( )
    ierr=1
    do while (ierr /= 0)
      read (unit, nml=segments_nml, iostat=io, end=40)
      ierr = check_nml_error (io, 'segments_nml')
    enddo
40  continue
    call close_file (unit)

  endif

!-----------------------------------------------------------------------
! write version number and namelists to logfile
!-----------------------------------------------------------------------

  call write_version_number (version, tag)
  if (mpp_pe() == mpp_root_pe()) write (stdlog(), nml=segments_nml)

!-----------------------------------------------------------------------
! write out overall beginning and ending times
!-----------------------------------------------------------------------

  call get_times ( Time_init, Time_end )
  call write_timestamps ( Time_init, Time_end, num_segs )

!-----------------------------------------------------------------------
! write out beginning and ending times for each segment
!-----------------------------------------------------------------------

!!$  if ( num_segs > 1 ) then
!!$     do nseg=1,num_segs
!!$        call get_times ( Time_init, Time_end, seg=nseg, nsegs=num_segs )
!!$        call write_timestamps ( Time_init, Time_end )
!!$     enddo
!!$  endif

!-----------------------------------------------------------------------
! write out lat/lon grid for grid_spec file creation
!-----------------------------------------------------------------------
  
!!$  call namelist_read (                                                 &
!!$                                                           ix, iy, iz, &
!!$                                                         nsteps, nbuf, &
!!$                             lxopen, lyopen, lcartesian, lcylindrical, &
!!$                                   rlonmin, rlonmax, rlatmin, rlatmax, &
!!$                                                           pbot, ptop, &
!!$                                                                 zeta )
!!$
!!$  allocate (zetaw(iz))
!!$  zetaw = zeta(1:iz)
!!$
!!$  unit = open_file ('horiz_grid.out', action='write')
!!$
!!$  if (mpp_pe() == mpp_root_pe()) write (unit,'(2i6,4f14.8)')           &
!!$                             ix, iy, rlonmin, rlonmax, rlatmin, rlatmax
  
  call close_file (unit)

!-----------------------------------------------------------------------
! exit mpi
!-----------------------------------------------------------------------

  call mpp_exit()
  
end program zetac_write_stuff

