module zetac_tracer_mod

!-----------------------------------------------------------------------
! module information
!-----------------------------------------------------------------------

use fms_mod,                     only : stdlog, stdout,                &
                                        write_version_number,          &
                                        mpp_pe, mpp_root_pe
use field_manager_mod,           only : model_atmos
use tracer_manager_mod,          only : register_tracers,              &
                                        get_tracer_names,              &
                                        get_tracer_index,              &
                                        query_method, NO_TRACER
use vert_advection_mod,          only : SECOND_CENTERED,               &
                                        FOURTH_CENTERED,               &
                                        VAN_LEER_LINEAR,               &
                                        FINITE_VOLUME_LINEAR,          &
                                        FINITE_VOLUME_PARABOLIC,       &
                                        FINITE_VOLUME_PARABOLIC2,      &
                                        FLUX_FORM,                     &
                                        ADVECTIVE_FORM

!-----------------------------------------------------------------------
! public interfaces
!-----------------------------------------------------------------------

implicit none
private
public get_info, get_method, tracer_type

character(len=*), parameter :: module='zetac_tracer_mod'
logical :: do_init=.true.

type tracer_type
  real, pointer, dimension(:,:,:,:) :: data  => NULL()
  real, pointer, dimension(:,:,:)   :: tend  => NULL()
  real, pointer, dimension(:,:,:)   :: fall  => NULL()
  real, pointer, dimension(:,:)     :: prec  => NULL()
  integer, pointer                  :: index => NULL()
end type tracer_type

!-----------------------------------------------------------------------
! version control information
!-----------------------------------------------------------------------

character(len=128)  :: version =  '$Id: zetac_tracer.f90,v 1.1.2.2.2.2 2004/10/19 22:58:19 stg Exp $'
character(len=128)  :: tag     =  '$Name:  $'

contains

!#######################################################################

subroutine get_info ( name, longname, units )

character(len=*), intent(in)  :: name
character(len=*), intent(out) :: longname, units

character(len=16) :: vname
integer :: index

  if (do_init) then 
     call tracer_init
  endif

  index = get_tracer_index ( model_atmos, name )
  if ( index == NO_TRACER ) then
     units = 'error'
     return
  endif

  call get_tracer_names ( model_atmos, index, vname, longname, units )

  return
end subroutine get_info

!#######################################################################

subroutine get_method ( index, method, name, num )

integer,                                          intent (in)   ::     &
                                                                index

character(len=*),                                 intent (in)   ::     &
                                                               method

character(len=*), optional,                       intent (out)  ::     &
                                                                 name

integer, optional,                                intent (out)  ::     &
                                                                  num

character(len=32) :: name_scheme
logical :: ldum

!-----------------------------------------------------------------------
! initialize module
!-----------------------------------------------------------------------

  if ( do_init ) then
     call tracer_init
  endif
  
!-----------------------------------------------------------------------
! query method associated with tracer
!-----------------------------------------------------------------------

  ldum = query_method ( method, model_atmos, index, name_scheme )

!-----------------------------------------------------------------------
! associate name_scheme with number
!-----------------------------------------------------------------------

  if ( PRESENT(name) ) then
     name = name_scheme
  endif

  if ( PRESENT(num) ) then
     if (method(1:6) == "advect") then
        if ( name_scheme(1:6) == 'fourth' ) then
           num = FOURTH_CENTERED
        else if ( name_scheme(1:8) == 'van_leer' ) then
           num = VAN_LEER_LINEAR
        else if ( name_scheme(1:18) == 'finite_volume_line' ) then
           num = FINITE_VOLUME_LINEAR
        else if ( name_scheme(1:18) == 'finite_volume_para' ) then 
           num = FINITE_VOLUME_PARABOLIC2  ! uses Lin limiters
        else
           num = SECOND_CENTERED
        endif
     endif
    
     if (method(1:8) == "equation") then
        if ( name_scheme(1:4) == 'flux' ) then
           num = FLUX_FORM
        else
           num = ADVECTIVE_FORM
        endif
     endif

  endif
  
  return
end subroutine get_method

!#######################################################################

subroutine tracer_init

integer :: nprog, ndiag, nfam, nvars

!-----------------------------------------------------------------------
! write version number to logfile
!-----------------------------------------------------------------------

  call write_version_number (version, tag)

!-----------------------------------------------------------------------
! allocate multiple tracer array and store tracer indices
!-----------------------------------------------------------------------

  call register_tracers ( model_atmos, nvars, nprog, ndiag, nfam )

  if (mpp_pe() == mpp_root_pe()) then
!!$     write (stdlog(), '("Number of family tracers =",i3)') nfam
!!$     write (stdlog(), '("Number of prognostic tracers =",i3)') nprog
!!$     write (stdlog(), '("Number of diagnostic tracers =",i3)') ndiag
     write (stdout(), '("Number of family tracers =",i3)') nfam
     write (stdout(), '("Number of prognostic tracers =",i3)') nprog
     write (stdout(), '("Number of diagnostic tracers =",i3)') ndiag
  endif

  do_init = .false.

  return
end subroutine tracer_init

!#######################################################################

end module zetac_tracer_mod
