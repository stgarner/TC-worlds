module zetac_time_filter_mod

!-----------------------------------------------------------------------
! module information
!-----------------------------------------------------------------------

use mpp_mod,            only : mpp_pe, mpp_root_pe, input_nml_file
use fms_mod,            only : file_exist,                             &
                               close_file, open_namelist_file,         &
                               write_version_number, check_nml_error,  &
                               stdlog, stdout

use zetac_horiz_grid_type_mod,   only : horiz_grid_type
use zetac_time_pointers_mod,     only : ntime

!-----------------------------------------------------------------------
! public interfaces
!-----------------------------------------------------------------------

implicit none
private
public time_filter, time_filter_init

! namelist
real    :: robfac=0.05

real    :: robfac1
logical :: do_forward

character(len=*), parameter :: module='zetac_time_filter_mod'

!-----------------------------------------------------------------------
! version control information
!-----------------------------------------------------------------------

character(len=128)  :: version =  '$Id: zetac_time_filter.f90,v 1.1.2.1.2.5 2005/08/06 23:57:34 stg Exp $'
character(len=128)  :: tag     =  '$Name:  $'

contains

!#######################################################################

subroutine time_filter (                                        Mgrid, &
                                                       uu, vv, ta, qv, &
                                   filt_uu, filt_vv, filt_ta, filt_qv, &
                                                past, present, future )

type (horiz_grid_type),                      intent (in)    ::         &
                                                                Mgrid

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked, ntime), intent (in)    ::         &
                                                       uu, vv, ta, qv

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),        intent (out)   ::         &
                                   filt_uu, filt_vv, filt_ta, filt_qv

integer,                                     intent (in)    ::         &
                                                past, present, future

!-----------------------------------------------------------------------
! filter variables at middle time step
!-----------------------------------------------------------------------

  call time_filter_single ( Mgrid, uu, filt_uu, past, present, future )
  call time_filter_single ( Mgrid, vv, filt_vv, past, present, future )
  call time_filter_single ( Mgrid, ta, filt_ta, past, present, future )
  call time_filter_single ( Mgrid, qv, filt_qv, past, present, future )

  return
end subroutine time_filter

!#######################################################################

subroutine time_filter_single (                                        &
                                                                Mgrid, &
                                                          field, filt, &
                                                past, present, future )

type (horiz_grid_type),                      intent (in)    ::         &
                                                                Mgrid

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked, ntime), intent (in)    ::         &
                                                                field

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),        intent (out)   ::         &
                                                                 filt

integer,                                     intent (in)    ::         &
                                                past, present, future

  if ( do_forward ) return

!-----------------------------------------------------------------------
! filter variable at middle time step
!-----------------------------------------------------------------------

  filt = robfac*(field(:,:,:,past) + field(:,:,:,future)               &
                                - 2.*field(:,:,:,present))

  return
end subroutine time_filter_single

!#######################################################################

subroutine time_filter_init ( Mgrid, ndt_forward )

type (horiz_grid_type),                     intent (in)    ::          &
                                                                Mgrid

integer,                                    intent (out)   ::          &
                                                          ndt_forward

integer :: io, ierr

namelist /time_filter_nml/ robfac, ndt_forward

!-----------------------------------------------------------------------
! read namelist
!-----------------------------------------------------------------------

  read (input_nml_file, time_filter_nml, iostat=io)
  ierr = check_nml_error(io,'time_filter_nml')

!-----------------------------------------------------------------------
! write version number and namelist to logfile
!-----------------------------------------------------------------------

  call write_version_number (version, tag)
  if (mpp_pe() == mpp_root_pe())                                       &
                                  write (stdlog(), nml=time_filter_nml)

!-----------------------------------------------------------------------
! define time filter coefficient
!-----------------------------------------------------------------------

  robfac1 = 1.0 - 2.0*robfac

  do_forward = ( ndt_forward > 0 )

  return
end subroutine time_filter_init

!#######################################################################

end module zetac_time_filter_mod
