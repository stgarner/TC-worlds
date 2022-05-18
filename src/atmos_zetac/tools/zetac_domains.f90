module zetac_domains_mod

!-----------------------------------------------------------------------
! module information
!-----------------------------------------------------------------------

use mpp_mod,         only : mpp_pe, mpp_root_pe, mpp_npes, mpp_sum,    &
                            input_nml_file
use fms_mod,         only : file_exist, write_version_number,          &
                            close_file, open_namelist_file,            &
                            check_nml_error, set_domain,               &
                            stdlog, stdout
use mpp_domains_mod, only : mpp_domains_init, mpp_define_domains,      &
                            mpp_domains_set_stack_size,                &
                            mpp_get_data_domain,                       &
                            mpp_get_compute_domain,                    &
                            cyclic_global_domain, domain2D

use zetac_horiz_grid_type_mod, only : horiz_grid_type,                 &
                                      horiz_grid_type_define

!-----------------------------------------------------------------------
! public interfaces
!-----------------------------------------------------------------------

implicit none
private
public get_domains

character(len=*), parameter :: module='zetac_domains_mod'
logical :: do_init=.true.
integer :: halosize, npes

integer, dimension(2) :: layout

!-----------------------------------------------------------------------
! version control information
!-----------------------------------------------------------------------

character(len=128)  :: version =  '$Id: zetac_domains.f90,v 1.1.2.4.2.3 2004/08/02 20:10:46 stg Exp $'
character(len=128)  :: tag     =  '$Name: latest $'

contains

!#######################################################################

subroutine get_domains (                                               &
                                                     ix, iy, iz, nbuf, &
                                                                Mgrid, &
                                                lxopen, lyopen, bgrid )

!-----------------------------------------------------------------------
! arguments
!-----------------------------------------------------------------------

integer,                                       intent (in)  ::         &
                                                     ix, iy, iz, nbuf
 
type (horiz_grid_type),                        intent(out)  ::         &
                                                                Mgrid

logical,                                       intent (in)  ::        &
                                                       lxopen, lyopen

logical, optional,                             intent (in)  ::         &
                                                                bgrid

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

type (domain2D) :: Domain
type (domain2D) :: Domain_plusxy 
type (domain2D) :: Domain_u, Domain_v, Domain_nohalo

logical :: bgrid1

integer ::                                                             &
                                         ibc, iec, jbc, jec, kbc, kec, &
                                         ibd, ied, jbd, jed, kbd, ked, &
                                                   ibu, ieu, jbv, jev, &
                                                   ibv, iev, jbu, jeu, &
                                                   ibh, ieh, jbh, jeh, &
                                                   ibg, ieg, jbg, jeg

  if ( do_init ) then
     call domains_init
  endif

  ibg = 1  - nbuf
  ieg = ix + nbuf
  jbg = 1  - nbuf
  jeg = iy + nbuf

  if ( present(bgrid) ) then
     bgrid1 = bgrid
  else
     bgrid1 = .false.
  endif
  
!-----------------------------------------------------------------------
! set global domain limits
!-----------------------------------------------------------------------
  
  if ( bgrid1 ) then

     ibu = 1
     ibv = 1
     if ( lxopen .and. ix > 1 ) then
        ieu = ix-1
     else
        ieu = ix
     endif
     iev = ieu

     jbv = 1
     jbu = 1
     if ( lyopen .and. iy > 1 ) then
        jev = iy-1
     else
        jev = iy
     endif
     jeu = jev

  else

     if ( lxopen .and. ix > 1 ) then
        ibu = 0
     else
        ibu = 1
     endif
     ibv = 1
     ieu = ix
     iev = ix

     if ( lyopen .and. iy > 1 ) then
        jbv = 0
     else
        jbv = 1
     endif
     jbu = 1
     jev = iy
     jeu = iy

  endif

!-----------------------------------------------------------------------
! set stack size
!-----------------------------------------------------------------------

  npes = mpp_npes()
  halosize = 4*ix*iy*iz/npes
  call mpp_domains_set_stack_size ( halosize )

!-----------------------------------------------------------------------
! get decomposition of special or extended grids
!-----------------------------------------------------------------------

  call get_decomp (                                                    &
                                                                 nbuf, &
                                                   ibg, ieg, jbg, jeg, &
                                                     .false., .false., &
                                                        Domain_plusxy )

  call get_decomp (                                                    &
                                                                    0, &
                                                   ibu, ieu, jbu, jeu, &
                                                     .false., .false., &
                                                             Domain_u )

  call get_decomp (                                                    &
                                                                    0, &
                                                   ibv, iev, jbv, jev, &
                                                     .false., .false., &
                                                             Domain_v )

  call get_decomp (                                                    &
                                                                    0, &
                                                         1, ix, 1, iy, &
                                                     .false., .false., &
                                                        Domain_nohalo )

  call set_domain (Domain_nohalo)  !physics decomposition

!-----------------------------------------------------------------------
! get decomposition of dyanamics grid
!-----------------------------------------------------------------------

  call get_decomp (                                                    &
                                                                 nbuf, &
                                                         1, ix, 1, iy, &
                                                       lxopen, lyopen, &
                                                               Domain )

  call mpp_get_data_domain ( Domain, ibd, ied, jbd, jed )

  call mpp_get_compute_domain ( Domain, ibc, iec, jbc, jec )

!-----------------------------------------------------------------------
! add on external halos for extrapolations at open boundaries
!-----------------------------------------------------------------------

  if ( lxopen .and. ibc == 1 ) then
     ibh = ibd
  else
     ibh = ibc-1 !stg
     ibu = ibc
  endif
  if ( lxopen .and. iec == ix ) then
     ieh = ied
     ieu = max(1, ix-1)
  else
     ieh = iec+1 !stg
     ieu = iec
  endif
  if ( lyopen .and. jbc == 1 ) then
     jbh = jbd
  else
     jbh = jbc-1 !stg
     jbv = jbc
  endif
  if ( lyopen .and. jec == iy ) then
     jeh = jed
     jev = max(1, iy-1)
  else
     jeh = jec+1 !stg
     jev = jec
  endif

!-----------------------------------------------------------------------
! define vertical grid indices
!-----------------------------------------------------------------------

  kbd = 1
  ked = iz
  kbc = kbd+1
  kec = ked

!-----------------------------------------------------------------------
! define derived-type variable 'Mgrid'
!-----------------------------------------------------------------------

  call horiz_grid_type_define (                                        &
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

  return
end subroutine get_domains

!#######################################################################

subroutine get_decomp (                                                &
                                                                 nbuf, &
                                         iwest, ieast, jsouth, jnorth, &
                                                       lxopen, lyopen, &
                                                               Domain )
 
integer,                                        intent (in)  ::        &
                                                                  nbuf
 
integer,                                         intent (in)  ::       &
                                          iwest, ieast, jsouth, jnorth

logical,                                        intent (in)  ::        &
                                                        lxopen, lyopen

type (domain2D),                                intent(out)  ::        &
                                                                Domain

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

integer :: xbegin, xend, ybegin, yend
integer, dimension(4) :: global

!-----------------------------------------------------------------------
! define domains
!-----------------------------------------------------------------------

  global = (/iwest,ieast,jsouth,jnorth/)

  if (lxopen .and. lyopen) then
    call mpp_define_domains( global, layout, Domain,                   &
        xhalo=nbuf, yhalo=nbuf )

  else if (lxopen .and. (.not.lyopen)) then
    call mpp_define_domains( global, layout, Domain,                   &
        yflags=cyclic_global_domain,                                   &
        xhalo=nbuf, yhalo=nbuf )

  else if ((.not.lxopen) .and. lyopen) then
    call mpp_define_domains( global, layout, Domain,                   &
        xflags=cyclic_global_domain,                                   &
        xhalo=nbuf, yhalo=nbuf )

  else
    call mpp_define_domains( global, layout, Domain,                   &
        xflags=cyclic_global_domain, yflags=cyclic_global_domain,      &
        xhalo=nbuf, yhalo=nbuf )
  endif

  return
end subroutine get_decomp

!#######################################################################

subroutine domains_init

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

integer :: io, unit, ierr

namelist /zetac_layout_nml/ layout

!-----------------------------------------------------------------------
! read namelist
!-----------------------------------------------------------------------

  read (input_nml_file, nml=zetac_layout_nml, iostat=io)
  ierr = check_nml_error(io,'zetac_layout_nml')

!-----------------------------------------------------------------------
! write version number and namelist to logfile
!-----------------------------------------------------------------------

  call write_version_number (version, tag)
  if (mpp_pe() == mpp_root_pe())                                       &
                                 write (stdlog(), nml=zetac_layout_nml)

  call mpp_domains_init

  do_init = .false.

  return
end subroutine domains_init

!#######################################################################

end module zetac_domains_mod
