module zetac_horiz_grid_type_mod

!-----------------------------------------------------------------------
! module information
!-----------------------------------------------------------------------

use fms_mod,         only : write_version_number

use mpp_domains_mod, only : domain2D

!-----------------------------------------------------------------------
! public interfaces
!-----------------------------------------------------------------------

implicit none
private
public horiz_grid_type, horiz_grid_type_define

type horiz_grid_type

  integer :: ix, iy, iz, nbuf
  logical :: lxopen, lyopen
  integer :: ibc, iec, jbc, jec, kbc, kec
  integer :: ibd, ied, jbd, jed, kbd, ked
  integer :: ibu, ieu, jbv, jev
  integer :: ibh, ieh, jbh, jeh
  integer :: ibg, ieg, jbg, jeg

  type (domain2D) :: Domain
  type (domain2D) :: Domain_plusxy
  type (domain2D) :: Domain_u, Domain_v, Domain_nohalo

end type horiz_grid_type

character(len=*), parameter :: module='zetac_horiz_grid_type_mod'
logical :: do_init=.true.

!-----------------------------------------------------------------------
! version control information
!-----------------------------------------------------------------------

character(len=128)  :: version =  '$Id: zetac_horiz_grid_type.f90,v 1.1.2.4 2004/05/21 15:46:02 ck Exp $'
character(len=128)  :: tag     =  '$Name: latest $'
 
contains

!#######################################################################

subroutine horiz_grid_type_define (                                    &
                                                     ix, iy, iz, nbuf, &
                                                       lxopen, lyopen, &
                                         ibc, iec, jbc, jec, kbc, kec, &
                                         ibd, ied, jbd, jed, kbd, ked, &
                                                   ibu, ieu, jbv, jev, &
                                                   ibh, ieh, jbh, jeh, &
                                                   ibg, ieg, jbg, jeg, &
                                                Domain, Domain_plusxy, &
                                    Domain_u, Domain_v, Domain_nohalo, &
                                                                Mgrid )
   
integer,                                            intent (in) ::     &
                                                     ix, iy, iz, nbuf
    
logical,                                            intent (in) ::     &
                                                       lxopen, lyopen

integer,                                            intent (in) ::     &
                                         ibc, iec, jbc, jec, kbc, kec, &
                                         ibd, ied, jbd, jed, kbd, ked, &
                                                   ibu, ieu, jbv, jev, &
                                                   ibh, ieh, jbh, jeh, &
                                                   ibg, ieg, jbg, jeg
 
type (domain2D),                                    intent (in)  ::    &
                                                Domain, Domain_plusxy, &
                                    Domain_u, Domain_v, Domain_nohalo

type (horiz_grid_type),                             intent (out) ::    &
                                                                Mgrid

!-----------------------------------------------------------------------
! initialize module
!-----------------------------------------------------------------------

  if (do_init) then
    call horiz_grid_type_init
  endif

  Mgrid%ix   = ix
  Mgrid%iy   = iy
  Mgrid%iz   = iz
  
  Mgrid%nbuf = nbuf

  Mgrid%lxopen = lxopen
  Mgrid%lyopen = lyopen

  Mgrid%ibc = ibc
  Mgrid%iec = iec
  Mgrid%jbc = jbc
  Mgrid%jec = jec
  Mgrid%kbc = kbc
  Mgrid%kec = kec

  Mgrid%ibd = ibd
  Mgrid%ied = ied
  Mgrid%jbd = jbd
  Mgrid%jed = jed
  Mgrid%kbd = kbd
  Mgrid%ked = ked

  Mgrid%ibu = ibu
  Mgrid%ieu = ieu
  Mgrid%jbv = jbv
  Mgrid%jev = jev

  Mgrid%ibh = ibh
  Mgrid%ieh = ieh
  Mgrid%jbh = jbh
  Mgrid%jeh = jeh

  Mgrid%ibg = ibg
  Mgrid%ieg = ieg
  Mgrid%jbg = jbg
  Mgrid%jeg = jeg
   
  Mgrid%Domain        = Domain
  Mgrid%Domain_plusxy = Domain_plusxy
  Mgrid%Domain_u      = Domain_u
  Mgrid%Domain_v      = Domain_v
  Mgrid%Domain_nohalo = Domain_nohalo

  return
end subroutine horiz_grid_type_define

!#######################################################################

subroutine horiz_grid_type_init

!-----------------------------------------------------------------------
! write version number to logfile
!-----------------------------------------------------------------------

  call write_version_number (version, tag)

  do_init = .false.
  
  return

end subroutine horiz_grid_type_init

!#######################################################################

end module zetac_horiz_grid_type_mod
