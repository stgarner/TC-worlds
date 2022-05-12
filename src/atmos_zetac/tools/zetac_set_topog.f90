module zetac_set_topog_mod

!-----------------------------------------------------------------------
! module information
!-----------------------------------------------------------------------

use netcdf 
use mpp_mod,        only : mpp_pe, mpp_root_pe, mpp_broadcast,         &
                           mpp_sum, input_nml_file
use mpp_io_mod,     only : MPP_NETCDF, MPP_MULTI,                      &
                           MPP_SINGLE, MPP_RDONLY,                     &
                           mpp_open, mpp_close, mpp_read,              &
                           axistype, fieldtype,                        &
                           mpp_get_atts, mpp_get_fields,               &
                           mpp_get_axes, mpp_get_info,                 &
                           mpp_get_axis_data
use fms_mod,        only : file_exist, close_file, open_namelist_file, &
                           write_version_number, mpp_pe, mpp_root_pe,  &
                           check_nml_error, stdlog, stdout,            &
                           error_mesg, FATAL, WARNING

use mpp_domains_mod,  only : mpp_get_data_domain
use horiz_interp_mod, only : horiz_interp_new, horiz_interp,           &
                             horiz_interp_type
use constants_mod,    only : radius, pi, radian

use zetac_domains_mod,           only : get_domains
use zetac_horiz_grid_type_mod,   only : horiz_grid_type
use zetac_horiz_metric_type_mod, only : horiz_metric_type
use zetac_update_halos_mod,      only : update_halos

!-----------------------------------------------------------------------
! public interfaces
!-----------------------------------------------------------------------

implicit none
private
public set_topog, nc_read_util

character(len=*), parameter :: module='zetac_set_topog_mod'
character (len=16) :: name
logical :: do_init=.true.

integer :: threading_rd=MPP_MULTI, fileset_rd=MPP_SINGLE

real, parameter :: tiny=1.0e-10

character*30, parameter :: topog_file='INPUT/topography_data.nc'
character*30, parameter :: sst_file='INPUT/sst_data.nc'

real, dimension(:), allocatable :: vlon
real, dimension(:), allocatable :: ulat

real :: dlon, dlat, dx, dy
integer :: ngridx, ngridy

!namelist:
logical :: read_topog=.true.
logical :: read_sst=.false.
real    :: grids=1.5, hclip=0.0

!-----------------------------------------------------------------------
! version control information
!-----------------------------------------------------------------------

character(len=128)  :: version =  '$Id: zetac_set_topog.f90,v 1.1.2.6.2.12 2005/06/16 18:22:01 stg Exp $'
character(len=128)  :: tag     =  '$Name:  $'

contains

!#######################################################################

subroutine set_topog (                                          Mgrid, &
                                                  Hmetric, topog, sst )

type (horiz_grid_type),                             intent (in)   ::   &
                                                                Mgrid

type (horiz_metric_type),                           intent (in)   ::   &
                                                              Hmetric

real, dimension(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed),              &
                                                    intent (out)  ::   &
                                                           topog, sst

!------------------------------------------------------------------------
! local allocations
!------------------------------------------------------------------------

real,dimension(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed) ::             &
                                                       scalex, scaley

real    :: scale
integer :: i, j
integer :: unit, io, ierr

namelist /topog_data_nml/ read_topog, read_sst, grids, hclip

!-----------------------------------------------------------------------
! read namelist
!-----------------------------------------------------------------------

  read (input_nml_file, nml=topog_data_nml, iostat=io)
  ierr = check_nml_error(io,'topog_data_nml')

!-----------------------------------------------------------------------
! write version number and namelist to logfile
!-----------------------------------------------------------------------
 
  if ( do_init ) call set_topog_init

  if (mpp_pe() == mpp_root_pe())                                       &
                                   write (stdlog(), nml=topog_data_nml)

!-----------------------------------------------------------------------
! grid geometry
!-----------------------------------------------------------------------

  dlon = Hmetric%dlon
  dlat = Hmetric%dlat
  dx   = Hmetric%dx
  dy   = Hmetric%dy

  allocate (ulat(Mgrid%jbd:Mgrid%jed))
  allocate (vlon(Mgrid%ibd:Mgrid%ied))
  ulat = Hmetric%ulat(Mgrid%jbd:Mgrid%jed)
  vlon = Hmetric%vlon(Mgrid%ibd:Mgrid%ied)

  if (vlon(Mgrid%ibh) < 0.0) vlon = vlon + 360.0

!-----------------------------------------------------------------------
! halo size for interpolations
!-----------------------------------------------------------------------

  scale = grids*sqrt(dx*dy)/radius

  do j=Mgrid%jbd,Mgrid%jed
     do i=Mgrid%ibd,Mgrid%ied
        scalex(i,j) = scale/cos(ulat(j)/radian)
        scaley(i,j) = scale
     enddo
  enddo

  ngridx = nint(maxval(scalex)*radian/dlon)
  ngridy = nint(maxval(scaley)*radian/dlat)

  ierr = 0.5*(Mgrid%ix - (Mgrid%iec - Mgrid%ibc + 1)) - 1.
  ngridx = min(ngridx, ierr)

!-----------------------------------------------------------------------
! get real or idealized topography
!-----------------------------------------------------------------------

  if ( read_topog ) then
     call interp_topog ( Mgrid, Hmetric, scalex, scaley, topog )
  else
     call create_topog ( Mgrid, Hmetric, topog )
  endif

  if ( read_sst ) then
     call input_sst ( Mgrid, Hmetric, sst )
  else
     call create_sst ( Mgrid, Hmetric, sst )
  endif

  call update_halos ( Mgrid, topog )
  call update_halos ( Mgrid, sst )

  if (hclip /= 0.0) topog = min(topog, hclip)

  return
end subroutine set_topog

!#######################################################################

subroutine interp_topog ( Mgrid, Hmetric, scalex, scaley, topog )

type (horiz_grid_type),                             intent (in)   ::   &
                                                                Mgrid

type (horiz_metric_type),                           intent (in)   ::   &
                                                              Hmetric

real, dimension(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed),              &
                                                    intent (in)   ::   &
                                                       scalex, scaley
   
real, dimension(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed),              &
                                                    intent (out)  ::   &
                                                                topog

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

real, dimension(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed) :: topo

integer, parameter :: nres=30
integer, parameter :: lx=nres*360,ly=nres*180  ! dimensions of topo file

real :: flonmin, flatmin, flonmax, flatmax
real :: dflon, dflat

type(fieldtype), allocatable :: Fields(:)

real, allocatable, dimension(:,:) :: htopog
real, allocatable, dimension(:)   :: flon, flat

logical :: fold

integer :: iounit, length, namelen
integer :: ib2, ie2, jb2, je2
integer :: ib1, ie1, jb1, je1
integer :: i, j
integer :: n, ndim, nvar, natt, ntimes, ier

  if (Mgrid%lxopen) then
     ib1 = max(Mgrid%ibg, Mgrid%ibd - ngridx)
     ie1 = min(Mgrid%ieg, Mgrid%ied + ngridx)
  else
     ib1 = mod(Mgrid%ibc - ngridx + 2*Mgrid%ix - 1, Mgrid%ix) + 1
     ie1 = mod(Mgrid%iec + ngridx + 2*Mgrid%ix - 1, Mgrid%ix) + 1
  endif  

  jb1 = max(Mgrid%jbg,   Mgrid%jbd - ngridy)
  je1 = min(Mgrid%jeg-1, Mgrid%jed + ngridy)

!-----------------------------------------------------------------------
! coordinates arrays for input topography (grid dimensions hard-wired)
!-----------------------------------------------------------------------

  allocate (flon(lx))
  allocate (flat(ly))

  flonmin =   0.0
  flonmax = 360.0
  flatmin = -90.0
  flatmax =  90.0

  dflon = (flonmax - flonmin)/lx
  dflat = (flatmax - flatmin)/ly

  do i=1,lx
    flon(i) = flonmin + (i - 0.5)*dflon
  enddo
  do j=1,ly
    flat(j) = flatmin + (j - 0.5)*dflat
  enddo

  if (Hmetric%ulat(jb1) < flat(1)) then
     do j=jb1+1,je1
        if (Hmetric%ulat(j) >= flat(1)) exit
     enddo
     jb1 = j
  endif
  if (Hmetric%ulat(je1) > flat(ly)) then
     do j=je1-1,jb1,-1
        if (Hmetric%ulat(j) <= flat(ly)) exit
     enddo
     je1 = j
  endif

!-----------------------------------------------------------------------
! map domains to external grid
!-----------------------------------------------------------------------

  call map_horiz_grid (                                                &
                                                  Hmetric, flon, flat, &
                               ib1, ie1, jb1, je1, ib2, ie2, jb2, je2, &
                                                                 fold )

!-----------------------------------------------------------------------
! read in topography
!-----------------------------------------------------------------------

  allocate (htopog(ib2:ie2,jb2:je2))

  if (fold) then
     n = lx
  else
     n = 0
  endif

  write (stdout(), '("zetac_set_topog: reading from ",a40)') topog_file

  namelen = len(trim(topog_file))
  call mpp_open (                     iounit, topog_file(1:namelen-3), &
                                   action=MPP_RDONLY, form=MPP_NETCDF, &
                             threading=MPP_SINGLE, fileset=MPP_SINGLE )

  call mpp_get_info ( iounit, ndim, nvar, natt, ntimes )

  allocate ( Fields(nvar) )
  call mpp_get_fields ( iounit, Fields )
!  do i=1,nvar
!     call mpp_get_atts ( Fields(i), name )
!  enddo
!stg  call mpp_close (iounit)

  call nc_read_util (                                                  &
                                          topog_file, 'topog', htopog, &
                                                n, ib2, ie2, jb2, je2 )
!htopog = max(0.0,htopog)  !stg

!-----------------------------------------------------------------------
! interpolate topography
!-----------------------------------------------------------------------

  call interp_xy (                                                     & 
                                            Mgrid, ib2, ie2, jb2, je2, &
                                               flon, flat, vlon, ulat, &
                                                       scalex, scaley, &
                                                         htopog, topo )

  deallocate ( flon, flat )
  deallocate ( htopog )

  do j=Mgrid%jbh,Mgrid%jeh
     do i=Mgrid%ibh,Mgrid%ieh
        topog(i,j) = max(0.0, topo(i,j))
     enddo
  enddo

  return
end subroutine interp_topog

!#######################################################################

subroutine input_sst ( Mgrid, Hmetric, sst )

type (horiz_grid_type),                             intent (in)   ::   &
                                                                Mgrid

type (horiz_metric_type),                           intent (in)   ::   &
                                                              Hmetric

real, dimension(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed),              &
                                                    intent (out)  ::   &
                                                                  sst

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

type(fieldtype), allocatable :: Fields_sst(:)
type(axistype),  allocatable :: Axes_sst(:)

integer :: n, ndim, nvar, natt, ntimes
integer :: unit_sst, length, namelen
integer :: ix, iy
logical :: found_sst

!-----------------------------------------------------------------------
! open sst data file
!-----------------------------------------------------------------------

  write (stdout(), '("zetac_set_topog: reading from ",a40)') sst_file

  namelen = len(trim(sst_file))
  call mpp_open (                     unit_sst, sst_file(1:namelen-3), &
                                   action=MPP_RDONLY, form=MPP_NETCDF, &
                           threading=threading_rd, fileset=fileset_rd )

!-----------------------------------------------------------------------
! get number of: dimensions, variables, attributes, times
!-----------------------------------------------------------------------

  call mpp_get_info ( unit_sst, ndim, nvar, natt, ntimes )

  allocate (Fields_sst(nvar))
  call mpp_get_fields ( unit_sst, Fields_sst )
  allocate (Axes_sst(ndim))
  call mpp_get_axes ( unit_sst, Axes_sst )

!-----------------------------------------------------------------------
! read in axis metadata
!-----------------------------------------------------------------------

  do n=1,2
     call mpp_get_atts ( Axes_sst(n), name, len=length )
     if ( name(1:3) == 'lon' ) then
        ix = length
     else if ( name(1:3) == 'lat' ) then
        iy = length
     endif
  enddo

  if (ix /= Mgrid%ix .or. iy /= Mgrid%iy) then
     call error_handler ( "sst has different grid", FATAL )
  endif

  do n=1,nvar
     call mpp_get_atts ( Fields_sst(n), name )
     found_sst = (name == "SST")
     if (found_sst) then
        call mpp_read ( unit_sst, Fields_sst(1), sst )
        exit
     endif
  enddo

  if (.not. found_sst) then
     call error_handler ( "sst variable not fouund", FATAL )
  endif

!stg  call mpp_close (unit_sst) 
  deallocate (Fields_sst, Axes_sst)

  return
end subroutine input_sst

!#######################################################################

subroutine map_help (                                                  &
                                                    Hmetric, lon, lat, &
                                   ib, ie, jb, je, ib2, ie2, jb2, je2 )

type (horiz_metric_type),           intent (in)            ::  Hmetric
real, dimension(:,:),               intent (in)            :: lon, lat
integer,                            intent (in)  :: ib,  ie,  jb,  je
integer,                            intent (out) :: ib2, ie2, jb2, je2

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

real    :: wlon, elon, slat, nlat
integer :: i, j, idim, jdim, i1, i2, ib1, ie1

  wlon = Hmetric%vlon(ib)
  elon = Hmetric%vlon(ie)
  slat = Hmetric%ulat(jb)
  nlat = Hmetric%ulat(je)

  idim = size(lon,1) ; jdim = size(lon,2)

  do j=1,jdim
     do i=1,idim
        if ( lon(i,j) >= wlon ) exit
     enddo
     i1 = max(1,i-1)
     do i=idim,1,-1
        if ( lon(i,j) <= elon ) exit
     enddo
     i2 = min(idim,i+1)
     if ( maxval(lat(i1:i2,j)) >= slat ) exit
     ib1 = i1 ; ie1 = i2
  enddo
  jb2 = max(1,j-1)
  if ( j == jb ) then
     ib1 = i1 ; ie1 = i2
  endif

  do j=jdim,1,-1
     do i=1,idim
        if ( lon(i,j) >= wlon ) exit
     enddo
     i1 = max(1,i-1)
     do i=idim,1,-1
        if ( lon(i,j) <= elon ) exit
     enddo
     i2 = min(idim,i+1)
     if ( minval(lat(i1:i2,j)) <= nlat ) exit
     ib2 = i1 ; ie2 = i2
  enddo
  je2 = min(jdim,j+1)
  if ( j == je ) then
     ib2 = i1 ; ie2 = i2
  endif

  ib2 = min(ib1,ib2) ; ie2 = max(ie1,ie2)

  return
end subroutine map_help

!#######################################################################

subroutine map_horiz_grid (                                            & 
                                                    Hmetric, lon, lat, &
                                   ib, ie, jb, je, ib2, ie2, jb2, je2, &
                                                                 fold )

type (horiz_metric_type), intent (in)  ::  Hmetric
real, intent(in), dimension(:) :: lon, lat
integer, intent(in)  :: ib,  ie,  jb,  je
integer, intent(out) :: ib2, ie2, jb2, je2
logical, intent(out) :: fold

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

real :: slat, nlat, wlon, elon, wlon_ext, elon_ext
integer :: idim, jdim

  wlon = Hmetric%vlon(ib)
  elon = Hmetric%vlon(ie)
  slat = Hmetric%ulat(jb)
  nlat = Hmetric%ulat(je)

  idim = size(lon) ; jdim = size(lat)

  wlon_ext = lon(1)
  elon_ext = lon(idim) - 360.0
  wlon = mod(wlon - wlon_ext + 360.0, 360.0) + wlon_ext
  elon = mod(elon - elon_ext + 360.0, 360.0) + elon_ext

  call get_indices ( idim, ib2, ie2, lon, wlon, elon )
  call get_indices ( jdim, jb2, je2, lat, slat, nlat )

  fold = ( ib2 > ie2 )
  if ( fold ) ie2 = ie2 + idim

  return
end subroutine map_horiz_grid

!#######################################################################

subroutine get_indices (                                               &
                                                     ndim, nbeg, nend, &
                                 axes_src, axes_des_beg, axes_des_end )

integer,                    intent (in)  :: ndim
integer,                    intent (out) :: nbeg, nend
real, dimension(ndim),      intent (in)  :: axes_src
real,                       intent (in)  :: axes_des_beg, axes_des_end

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

  integer :: n

  nbeg = -999
  do n=1,ndim
     if ( axes_src(n) < axes_des_beg + tiny ) nbeg = n
  enddo

  nend = -999
  do n=ndim,1,-1
     if ( axes_src(n) > axes_des_end - tiny ) nend = n
  enddo

  if ( nbeg == -999 .or. nend  == -999 ) call error_handler            &
                  ('destination grid is not inside source grid', FATAL) 

  return
end subroutine get_indices

!#######################################################################

subroutine interp_xy (                                                 &
                                            Mgrid, lx1, lx2, ly1, ly2, &
                                               flon, flat, clon, clat, &
                                                       scalex, scaley, &
                                                         fdata, cdata )

type (horiz_grid_type),                          intent (in)  ::       &
                                                                Mgrid

integer,intent(in) :: lx1, lx2, ly1, ly2

real, dimension(:),                              intent(in)  ::  flon
real, dimension(:),                              intent(in)  ::  flat
real, dimension(Mgrid%ibd:Mgrid%ied),            intent(in)  ::  clon
real, dimension(Mgrid%jbd:Mgrid%jed),            intent(in)  ::  clat

real, dimension(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed),              &
                                                 intent(in)  ::        &
                                                       scalex, scaley

real, dimension(lx1:lx2,ly1:ly2),                intent(in)  ::        &
                                                                fdata

real, dimension(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed),              &
                                                 intent(out) ::        &
                                                                cdata

!------------------------------------------------------------------------
! local allocations
!------------------------------------------------------------------------

integer :: i, i1, ibh, ieh
integer :: j, j1, jbh, jeh

real,dimension(lx1:lx2)                     ::                    flam
real,dimension(ly1:ly2)                     ::  cosfphi, sinfphi, fphi
real,dimension(Mgrid%ibd:Mgrid%ied)         ::                    clam
real,dimension(Mgrid%jbd:Mgrid%jed)         ::  coscphi, sincphi, cphi
real,dimension(Mgrid%ibd:Mgrid%ied,lx1:lx2) ::                 cosdlam

real, dimension(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed) ::  suma,sumb

real    :: dfphi, dflam
real    :: cj, sj, cj1, sj1, ccj, ssj, arc, arc1, w, sum
integer :: idiff, jdiff, i0, j0, ia, ib, ja, jb
integer :: idim, lx

  ibh = Mgrid%ibh
  ieh = Mgrid%ieh
  jbh = Mgrid%jbh
  jeh = Mgrid%jeh

  idim = size(flon)

  lx = min(lx2,idim)
  do i=lx1,lx
     flam(i) = flon(i)/radian
  enddo
  if (lx2 > lx) then
     do i=1,lx2-lx
        flam(i+lx) = (flon(i) + 360.)/radian
     enddo
  endif
  do i=ibh,ieh
     clam(i) = clon(i)/radian
  enddo

  if (clam(ibh) < flam(lx1)) clam = clam + 2.*pi
  do j=ly1,ly2
     fphi(j) = flat(j)/radian
     cosfphi(j) = cos(fphi(j))
     sinfphi(j) = sin(fphi(j))
  enddo
  do j=jbh,jeh
     cphi(j) = clat(j)/radian
     coscphi(j) = cos(cphi(j))
     sincphi(j) = sin(cphi(j))
  enddo

  do i=ibh,ieh
     do i1=lx1,lx2
        cosdlam(i,i1) = cos(flam(i1) - clam(i))
     enddo
  enddo

  dflam = flam(lx1+1) - flam(lx1)
  dfphi = fphi(ly1+1) - fphi(ly1)

  suma = tiny
  sumb = tiny

  do j=jbh,jeh
     if (j < 1 .or. j > Mgrid%iy) cycle
     cj = coscphi(j)
     sj = sincphi(j)
     do i=ibh,ieh
        idiff = int( scalex(i,j)/dflam )
        jdiff = int( scaley(i,j)/dfphi )
        i0 = int( (clam(i) - flam(lx1))/dflam + tiny ) + lx1
        j0 = int( (cphi(j) - fphi(ly1))/dfphi + tiny ) + ly1
        ia = max( lx1, i0 - idiff )
        ib = min( lx2, i0 + idiff )
        ja = max( ly1, j0 - jdiff )
        jb = min( ly2, j0 + jdiff )
        do j1=ja+1,jb
           cj1 = cosfphi(j1)
           sj1 = sinfphi(j1)
           ccj = cj*cj1
           ssj = sj*sj1
           do i1=ia+1,ib
              arc = acos( min(1.0, ccj*cosdlam(i,i1) + ssj) )
              arc1 = arc/max(arc, scaley(i,j))
              w = cj1*(0.42 + 0.50*cos(pi*arc1) + 0.08*cos(2*pi*arc1))
              suma(i,j) = suma(i,j) + w*fdata(i1,j1)
              sumb(i,j) = sumb(i,j) + w
           enddo
        enddo
     enddo
  enddo

  do j=jbh,jeh
     do i=ibh,ieh
        cdata(i,j) = suma(i,j)/sumb(i,j)
     enddo
  enddo

  if (jbh < 1) then
  do j=jbh,0
     cdata(:,j) = cdata(:,1)
  enddo
  endif
  if (jeh > Mgrid%iy) then
  do j=Mgrid%iy+1,jeh
     cdata(:,j) = cdata(:,Mgrid%iy)
  enddo
  endif

  return
end subroutine interp_xy

!#######################################################################

subroutine create_topog ( Mgrid, Hmetric, topog )

type (horiz_grid_type),               intent (in)            ::        &
                                                                Mgrid

type (horiz_metric_type),             intent (in)            ::        &
                                                              Hmetric

real, dimension(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed),              &
                                      intent(out)            :: topog
				      
!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

real, dimension(Mgrid%jbh:Mgrid%jeh) :: dxu

integer :: unit, io, ierr
integer :: i, ibh, ieh
integer :: j, jbh, jeh

real :: bell, errf, cos2
real :: dx2, dy2, dxfac, dyfac
real :: alpha, csal, snal, fxy
real :: xp, xr, xr1, yp, yr, yr1

!namelist:
real :: hmnt=0.0e3, mntlon, mntlat
real :: xscale=0.0e3, yscale=0.0e3, xridge=0.0e3, yridge=0.0e3
real :: angle=0.0

namelist /topog_nml/ hmnt, mntlon, mntlat,                             &
                     xscale, yscale, xridge, yridge, angle
  
!-----------------------------------------------------------------------
! statement functions for idealized topography
!-----------------------------------------------------------------------

  bell (xp) = 1.0/(1.0 + xp*xp)
  errf (xp) = exp(-xp*xp)
  cos2 (xp) = 0.5 + 0.5*cos(pi*min(xp,1.0))

!-----------------------------------------------------------------------
! read namelist
!-----------------------------------------------------------------------

  read (input_nml_file, nml=topog_nml, iostat=io)
  ierr = check_nml_error(io,'topog_nml')

!-----------------------------------------------------------------------
! write version number and namelist to logfile
!-----------------------------------------------------------------------
 
  if (mpp_pe() == mpp_root_pe())                                       &
                                        write (stdlog(), nml=topog_nml)

  if ( hmnt == 0.0 .or. (xscale == 0.0 .and. yscale == 0.0) ) then
     topog = hmnt
     return
  endif

  ibh = Mgrid%ibh
  ieh = Mgrid%ieh
  jbh = Mgrid%jbh
  jeh = Mgrid%jeh

  dxu  = Hmetric%dxu(Mgrid%jbh:Mgrid%jeh)

!-----------------------------------------------------------------------
! define idealized topography: symmetric ridge with half-widths 'xscale' 
! and 'yscale', plateau half-length 'ridge', and orientation 'angle' 
! deg clockwise from north-south
!-----------------------------------------------------------------------
  
  dx2 = (xscale/radian*radius)**2
  dy2 = (yscale/radian*radius)**2
  alpha = angle/radian
  dxfac = dx2/(dx2 + tiny)
  dyfac = dy2/(dy2 + tiny)
  csal = cos(alpha) * dxfac
  snal = sin(alpha) * dyfac

  do j=jbh,jeh
     yp = ((Hmetric%rlatmin - mntlat)/dlat + j-0.5)*dy
     do i=ibh,ieh
        xp  = ((Hmetric%rlonmin - mntlon)/dlon + i-0.5)*dxu(j)
        xr  = csal*xp - snal*yp
        yr  = snal*xp + csal*yp
        xr1 = max( dim(xr,xridge), dim(-xridge,xr) )
        yr1 = max( dim(yr,yridge), dim(-yridge,yr) )
        fxy = sqrt(xr1**2/(dx2 + tiny) + yr1**2/(dy2 + tiny))
        topog(i,j) = hmnt*cos2(fxy)
     enddo
  enddo

!-----------------------------------------------------------------------
! in case of cyclic boundaries:
!-----------------------------------------------------------------------

  if ( .not. Mgrid%lxopen .and. xscale < radius ) then
     do i=ibh,ieh
        xp = max(0.0, min(1.0, float(i-1)/Mgrid%ix))
        xp = (0.5 - 0.5*cos(2.0*pi*xp))
        do j=jbh,jeh
           topog(i,j) = topog(i,j)*xp
        enddo
     enddo
  endif
  if ( .not. Mgrid%lyopen .and. yscale < radius ) then
     do j=jbh,jeh
        yp = max(0.0, min(1.0, float(j-1)/Mgrid%iy))
        yp = (0.5 - 0.5*cos(2.0*pi*yp))
        do i=ibh,ieh
           topog(i,j) = topog(i,j)*yp
        enddo
     enddo
  endif

  return
end subroutine create_topog

!#######################################################################

subroutine create_sst ( Mgrid, Hmetric, sst )

type (horiz_grid_type),               intent (in)            ::        &
                                                                Mgrid

type (horiz_metric_type),             intent (in)            ::        &
                                                              Hmetric

real, dimension(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed),              &
                                      intent(out)            ::   sst
				      
!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

real, dimension(Mgrid%jbh:Mgrid%jeh) :: dxu

integer :: unit, io, ierr
integer :: i, ibh, ieh
integer :: j, jbh, jeh

real :: dx2, dy2, dxfac, dyfac
real :: xp, xr, xr1, yp, yr, yr1

!namelist:
logical :: use_constant=.true.
real :: sst_const=300.
real :: xscale=0.0e3, yscale=0.0e3

namelist /sst_nml/ use_constant, sst_const,                            &
                   xscale, yscale
  
!-----------------------------------------------------------------------
! read namelist
!-----------------------------------------------------------------------

  read (input_nml_file, nml=sst_nml, iostat=io)
  ierr = check_nml_error(io,'sst_nml')

!-----------------------------------------------------------------------
! write version number and namelist to logfile
!-----------------------------------------------------------------------
 
  if (mpp_pe() == mpp_root_pe())                                       &
                                         write (stdlog(), nml=sst_nml)

  if ( use_constant ) then
     sst = sst_const
     return
  endif

  ibh = Mgrid%ibh
  ieh = Mgrid%ieh
  jbh = Mgrid%jbh
  jeh = Mgrid%jeh

  dxu  = Hmetric%dxu(Mgrid%jbh:Mgrid%jeh)

!-----------------------------------------------------------------------
! define idealized sst
!-----------------------------------------------------------------------
  
  dx2 = (xscale/radian*radius)**2
  dy2 = (yscale/radian*radius)**2
  dxfac = dx2/(dx2 + tiny)
  dyfac = dy2/(dy2 + tiny)

  do j=jbh,jeh
     yp = 0.
     do i=ibh,ieh
        xp  = 0.
        sst(i,j) = sst_const
     enddo
  enddo

  return
end subroutine create_sst

!#######################################################################

subroutine nc_read_util ( filename, name, var, lx, ib, ie, jb, je )

character(len=*),                                    intent (in)    :: &
                                                        filename, name

real, dimension(:,:),                                intent (out)   :: &
                                                                   var

integer,                                             intent (in)    :: &
                                                    lx, ib, ie, jb, je

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

real, allocatable, dimension(:,:) :: var_west, var_east

integer :: status, rstatus
integer :: NCID, VARID
integer :: start(3), count(3)
integer :: lxwest, lxeast

  status = NF90_OPEN (filename, NF90_NOWRITE, NCID)

  if ( status /= 0 ) then
     call error_handler ('netcdf file not found', FATAL) 
  endif

  status = NF90_INQ_VARID (NCID, name, VARID)

  if ( status /= 0 ) then
     call error_handler ('netcdf variable not found', FATAL) 
  endif

  if ( lx == 0 ) then

     start(1) = ib
     start(2) = jb
     start(3) = 1
     count(1) = size(var,1)
     count(2) = size(var,2)
     count(3) = 1

     rstatus = NF90_GET_VAR ( NCID, VARID, var, start, count )

     if ( rstatus /= 0 ) then
        call error_handler ('read failed', FATAL)
     endif

  else

     lxeast = lx - ib + 1
     lxwest = ie - lx

     allocate ( var_west(lxwest, size(var,2)) )
     allocate ( var_east(lxeast, size(var,2)) )

     start(1) = 1
     start(2) = jb
     start(3) = 1
     count(1) = lxwest
     count(2) = size(var,2)  
     count(3) = 1

     rstatus = NF90_GET_VAR ( NCID, VARID, var_west, start, count )

     if ( rstatus /= 0 ) then
        call error_handler ('read failed east of fold', FATAL) 
     endif

     var(lxeast+1:lxeast+lxwest,:) = var_west

     start(1) = ib
     count(1) = lxeast

     rstatus = NF90_GET_VAR ( NCID, VARID, var_east, start, count )

     if ( rstatus /= 0 ) then
        call error_handler ('read failed west of fold', FATAL)
     endif

     var(1:lxeast,:) = var_east

     deallocate ( var_west, var_east )

  endif

  status = NF90_CLOSE (NCID)

  return
end subroutine nc_read_util

!#######################################################################

subroutine nc_read_global ( filename, name, var )

character(len=*),                                    intent (in)    :: &
                                                        filename, name

real, dimension(:,:),                                intent (out)   :: &
                                                                   var

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

integer :: status, rstatus
integer :: NCID, VARID
integer :: start(3), count(3)

  status = NF90_OPEN (filename, NF90_NOWRITE, NCID)

  if ( status /= 0 ) then
     call error_handler ('netcdf file not found', FATAL) 
  endif

  status = NF90_INQ_VARID (NCID, name, VARID)

  if ( status /= 0 ) then
     call error_handler ('netcdf variable not found', FATAL)
  endif

  start(1) = 1
  start(2) = 1
  start(3) = 1
  count(1) = size(var,1)
  count(2) = size(var,2)
  count(3) = 1

  rstatus = NF90_GET_VAR ( NCID, VARID, var, start, count )

  if ( rstatus /= 0 ) then
     call error_handler ('netcdf global read failed', FATAL)
  endif

  status = NF90_CLOSE (NCID)

  return
end subroutine nc_read_global

!#######################################################################

subroutine set_topog_init
 
!-----------------------------------------------------------------------
! write namelist to logfile
!-----------------------------------------------------------------------

  call write_version_number (version, tag)

  do_init = .false.

  return
end subroutine set_topog_init

!#######################################################################

subroutine error_handler ( message, level )

character(len=*), intent(in) :: message
integer, intent(in) :: level

  call error_mesg ( module, message, level )

  return
end subroutine error_handler

!#######################################################################

end module zetac_set_topog_mod
