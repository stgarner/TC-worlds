module zetac_namelist_mod

use fms_mod,    only : write_version_number, mpp_pe, mpp_root_pe,      &
                       check_nml_error, stdlog, stdout
use mpp_mod,    only : input_nml_file

!-----------------------------------------------------------------------
! public interfaces
!-----------------------------------------------------------------------

implicit none
private
public namelist_read

!-----------------------------------------------------------------------
! version control information
!-----------------------------------------------------------------------

character(len=128)  :: version =  '$Id: zetac_namelist.f90,v 1.1.2.2.2.4 2004/09/29 16:09:17 stg Exp $'
character(len=128)  :: tag     =  '$Name:  $'

contains

!#######################################################################

subroutine namelist_read (                                             &
                                                           ix, iy, iz, &
                                                         nsteps, nbuf, &
                             lxopen, lyopen, lcartesian, lcylindrical, &
                                   rlonmin, rlonmax, rlatmin, rlatmax, &
                                                           pbot, ptop, &
                                                                 zeta )

integer,                                    intent (out)   ::          &
                                                           ix, iy, iz

integer,                                    intent (out)   ::          &
                                                         nsteps, nbuf

logical,                                    intent (out)   ::          &
                             lxopen, lyopen, lcartesian, lcylindrical

real,                                       intent (out)   ::          &
                                   rlonmin, rlonmax, rlatmin, rlatmax

real,                                       intent (out)   ::          &
                                                           pbot, ptop

real, dimension(20),                        intent (inout) ::          &
                                                                 zeta

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

integer :: io, ierr

namelist /grid_nml/       &
  ix,                     &
  iy,                     &
  iz,                     &
  lxopen,                 &
  lyopen,                 &
  lcartesian,             &
  lcylindrical,           &
  nsteps,                 &
  nbuf                       

namelist /region_nml/     &
  rlonmin,                &
  rlonmax,                &
  rlatmin,                &
  rlatmax

namelist /vertical_nml/   &
  pbot,                   &
  ptop,                   &
  zeta

!-----------------------------------------------------------------------
! read namelists
!-----------------------------------------------------------------------

  read (input_nml_file, nml=grid_nml, iostat=io)
  ierr = check_nml_error(io,'grid_nml')
  read (input_nml_file, nml=region_nml, iostat=io)
  ierr = check_nml_error(io,'region_nml')
  read (input_nml_file, nml=vertical_nml, iostat=io)
  ierr = check_nml_error(io,'vertical_nml')

!-----------------------------------------------------------------------
! write version number and namelists to logfile
!-----------------------------------------------------------------------

  call write_version_number (version, tag)

  if (mpp_pe() == mpp_root_pe())                                       &
                                         write (stdlog(), nml=grid_nml)
  if (mpp_pe() == mpp_root_pe())                                       &
                                       write (stdlog(), nml=region_nml)
  if (mpp_pe() == mpp_root_pe())                                       &
                                     write (stdlog(), nml=vertical_nml)

  return
end subroutine namelist_read

!#######################################################################

end module zetac_namelist_mod

