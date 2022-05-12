module zetac_history_mod

!-----------------------------------------------------------------------
! module information
!-----------------------------------------------------------------------

use mpp_mod, only : mpp_pe, mpp_root_pe, mpp_npes, mpp_sum,            &
                    mpp_clock_id, mpp_clock_begin,                     &
                    mpp_clock_end, input_nml_file
use fms_mod, only : file_exist, close_file, open_namelist_file,        &
                    write_version_number, check_nml_error,             &
                    stdlog, stdout, error_mesg, FATAL, WARNING
use constants_mod,               only : grav, cp_air, hlv
use diag_manager_mod,            only : register_static_field,         &
                                        register_diag_field,           &
                                        send_data, need_data
use time_manager_mod,            only : time_type, get_time,           &
                                        operator(+)

use zetac_axes_mod,              only : axes_history, axid,            &
                                        gridu, gridv, gridm, gridw
use zetac_axis_names_mod,        only : xu_, yu_, xv_, yv_, yg_,       &
                                        xm_, ym_, zm_, zw_
use zetac_coriolis_param_mod,    only : coriolis_param
use zetac_cgrid_interp_mod,      only : cgrid_interp_mass
use zetac_convert_var_mod,       only : get_pp, get_sd
use zetac_field_names_mod
use zetac_horiz_grid_type_mod,   only : horiz_grid_type
use zetac_horiz_metric_mod,      only : CYLINDRICAL
use zetac_horiz_metric_type_mod, only : horiz_metric_type
use zetac_ncdf_io_mod,           only : ncdf_fms_write
use zetac_phys_con_mod,          only : tref
use zetac_time_pointers_mod,     only : ntime
use zetac_tracer_mod,            only : get_info
use zetac_vert_metric_type_mod,  only : vert_metric_type

!-----------------------------------------------------------------------
! public interfaces
!-----------------------------------------------------------------------

implicit none
private
public write_history, write_history_init

integer :: zetac_write_history_clock
integer :: zetac_history_stdout_clock

integer, parameter :: num_surf_fields=3

integer :: id_vtot, id_um, id_vm, id_hdiv, id_vort, id_fz
integer :: id_ekin, id_epot, id_ekin_int, id_epot_int
integer :: id_kedt, id_pedt, id_kedt_int, id_pedt_int

integer, dimension(num_refs_fields) :: refs_id
integer, dimension(num_grid_fields) :: grid_id
integer, dimension(num_prog_fields) :: prog_id
integer, dimension(num_diag_fields) :: diag_id
integer, dimension(num_prog_fields) :: pert_id    ! use "prog" indices

integer, dimension(3) :: axid_uum, axid_vvm, axid_mmm, axid_mmw
integer, dimension(2) :: axid_mm, axid_gm

character (len=*),  parameter :: module='zetac_history_mod'
character (len=16), parameter :: model='atmos_mod'
character (len=8)  :: name
character (len=64) :: longname, units

integer            :: nsecs

real, parameter :: silly=-99999.9
real, parameter :: tiny=1.0e-10
logical, allocatable, dimension(:,:) :: p_is_min

real, allocatable, dimension(:,:)   :: udp, vdp, kdp, pdp, ktdp, ptdp
real, allocatable, dimension(:,:)   :: epref
real, allocatable, dimension(:,:,:) :: pp, up, tp, qp, gp
real, allocatable, dimension(:,:,:) :: ph, pf, sd
real, allocatable, dimension(:,:,:) :: vtot, hdiv, vort, ekin, epot
real, allocatable, dimension(:,:,:) :: ukin, vkin

real, allocatable, dimension(:)     :: zetam, zetaw, dxu
real, allocatable, dimension(:,:)   :: topog

integer, allocatable, dimension(:,:) :: loc
integer, allocatable, dimension(:)   :: offset
real,    allocatable, dimension(:)   :: var, tmp
real,    allocatable, dimension(:)   :: ulat, vlon, cosulat, cosvlat

real :: var1, var2, varbar, area, ainf, hinf=1.0e5
real :: pbot, ptop, delp, maxw, esave=0., qsave=0., detot, dqtot
real :: dy

logical :: first_time=.true.

!namelist:
integer :: nsteps_stdout=0
logical :: verbose=.true.
real    :: hmask=0.0, pcrit=940.e2
real    :: varsum=0., dqsum=0., qbar
integer :: nsum=0

integer :: ibc, iec, ibd, ied
integer :: jbc, jec, jbd, jed, iy
integer :: kbd, ked, iz

real, allocatable, dimension(:) :: temp, temp1, temp2

!-----------------------------------------------------------------------
! version control information
!-----------------------------------------------------------------------

character(len=128)  :: version =  '$Id: zetac_history.f90,v 1.1.2.10.2.9 2005/08/07 00:01:45 stg Exp $'
character(len=128)  :: tag     =  '$Name:  $'

contains

!#######################################################################

subroutine write_history (                                             &
                                      Mgrid, Hmetric, Time, Time_step, &
                        ps, dpu, dpv, dpm, uu, vv, oo, ta, qv, qc, sm, &
                                             om, um, vm, ttend, qtend, &
                                               rh, gz, fz, kedt, pedt, &
                                    ppref, uuref, taref, qvref, gzref, &
                                         elapsed_time, percent, mxval, &
                                                                 delt )

type (horiz_grid_type),                         intent (in)    ::      &
                                                                 Mgrid

type (horiz_metric_type),                       intent (in)    ::      &
                                                               Hmetric

type (time_type),                               intent (in)    ::      &
                                                       Time, Time_step

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed),             &
                                                intent (in)    ::      &
                                                                    ps

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),           intent (in)    ::      &
                             dpu, dpv, dpm, uu, vv, oo, ta, qv, qc, sm

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),           intent (inout) ::      &
                                                            om, um, vm

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),           intent (in)    ::      &
                                              ttend, qtend, rh, gz, fz

real, dimension(Mgrid%ibd:Mgrid%ied, Mgrid%jbd:Mgrid%jed,              &
                Mgrid%kbd:Mgrid%ked),           intent (inout) ::      &
                                                            kedt, pedt

real, dimension(Mgrid%jbg:Mgrid%jeg, Mgrid%iz), intent (in)    ::      &
                                     ppref, uuref, taref, qvref, gzref
				       
real,                                           intent (in)    ::      &
				          elapsed_time, percent, mxval

real,                                           intent (in)    :: delt

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

type (time_type)                                               ::      &
                                                            Time_next

logical :: used
real    :: dxj, dyj, cosj
real    :: minlon, minlat, maxlon, maxlat
integer :: i, j, k
integer :: n=0

!-----------------------------------------------------------------------
! write zee(x,y,z) and static fields f(z) and f(x,y)
!-----------------------------------------------------------------------

  call mpp_clock_begin ( zetac_write_history_clock )

  if ( first_time ) then

     call write_static ( Mgrid, Time, ppref, uuref, taref, qvref )

     do k=kbd+1,ked
        do j=jbd,jed
           epref(j,k) = taref(j,k)*(ppref(j,k) - ppref(j,k-1))
        enddo
     enddo
     epref(:,kbd) = epref(:,kbd+1)

     first_time = .false.

  endif

!-----------------------------------------------------------------------
! write prognostic fields f(x,y,z,t) 
!-----------------------------------------------------------------------

  call ncdf_fms_write ( Mgrid, Time, prog_id(uu_), gridu, uu )
  call ncdf_fms_write ( Mgrid, Time, prog_id(vv_), gridv, vv )
  call ncdf_fms_write ( Mgrid, Time, prog_id(oo_), gridw, oo )
  call ncdf_fms_write ( Mgrid, Time, prog_id(ta_), gridm, ta )
  call ncdf_fms_write ( Mgrid, Time, prog_id(qv_), gridm, qv )
  call ncdf_fms_write ( Mgrid, Time, prog_id(qc_), gridm, qc )

!-----------------------------------------------------------------------
! write diagnostic fields f(x,y,z)
!-----------------------------------------------------------------------

  Time_next = Time + Time_step

  call ncdf_fms_write ( Mgrid, Time, id_um, gridm, um )
  call ncdf_fms_write ( Mgrid, Time, id_vm, gridm, vm )

  vtot = sqrt(um**2 + vm**2)
  call ncdf_fms_write ( Mgrid, Time, id_vtot, gridm, vtot )

  ukin = 0.5*uu*uu*dpu
  vkin = 0.5*vv*vv*dpv
  call cgrid_interp_mass ( Mgrid, iz, ukin, vkin, um, vm )

  ekin = (um + vm) / dpm * (1. + qv)   ! kinetic energy
  epot = Cp_air*ta + Hlv*qv            ! internal energy   

  epot = epot - Cp_air*tref            ! adjust internal energy   

  do j=jbc,jec
     do i=ibc,iec
        kdp(i,j) = sum ( ekin(i,j,kbd+1:ked)*dpm(i,j,kbd+1:ked) )/Grav
        pdp(i,j) = sum ( epot(i,j,kbd+1:ked)*dpm(i,j,kbd+1:ked) )/Grav
        ktdp(i,j) = sum ( kedt(i,j,kbd+1:ked)*dpm(i,j,kbd+1:ked) )/Grav
        ptdp(i,j) = sum ( pedt(i,j,kbd+1:ked)*dpm(i,j,kbd+1:ked) )/Grav
     enddo
  enddo

  call ncdf_fms_write ( Mgrid, Time, id_ekin, gridm, ekin )
  call ncdf_fms_write ( Mgrid, Time, id_epot, gridm, epot )
  call ncdf_fms_write ( Mgrid, Time, id_ekin_int, gridm, kdp )
  call ncdf_fms_write ( Mgrid, Time, id_epot_int, gridm, pdp )

  call ncdf_fms_write ( Mgrid, Time, id_kedt, gridm, kedt )
  call ncdf_fms_write ( Mgrid, Time, id_pedt, gridm, pedt )
  call ncdf_fms_write ( Mgrid, Time, id_kedt_int, gridm, ktdp )
  call ncdf_fms_write ( Mgrid, Time, id_pedt_int, gridm, ptdp )

  if ( need_data(id_hdiv,Time_next) ) then
     do k=kbd,ked
        do j=jbc-1,jec
           cosj = cosvlat(j)
           do i=ibc-1,iec
              udp(i,j) = uu(i,j,k)*dpu(i,j,k)
              vdp(i,j) = vv(i,j,k)*dpv(i,j,k)*cosj
           enddo
        enddo
        do j=jbc,jec
           dxj = dxu(j)
           dyj = dy*cosulat(j)
           do i=ibc,iec
              hdiv(i,j,k) = ((udp(i,j) - udp(i-1,j))/dxj               &
                           + (vdp(i,j) - vdp(i,j-1))/dyj) / Grav
           enddo
        enddo
     enddo
     call ncdf_fms_write ( Mgrid, Time, id_hdiv, gridm, hdiv )
  endif

  if ( need_data(id_vort,Time_next) ) then
     do j=jbc,jec
        dxj = dxu(j)
        do i=ibc,iec
           vort(i,j,:) = 0.25*((vv(i+1,j,:) - vv(i-1,j,:)              &
                             + (vv(i+1,j-1,:) - vv(i-1,j-1,:)))/dxj    &
                            - ((uu(i,j+1,:) - uu(i,j-1,:))             &
                             + (uu(i-1,j+1,:) - uu(i-1,j-1,:)))/dy)
        enddo
     enddo
     call ncdf_fms_write ( Mgrid, Time, id_vort, gridm, vort )
  endif

  call get_pp ( ps, zetaw, ptop, ph )
  call get_pp ( ps, zetam, ptop, pf )

  call get_sd ( ta, pf, sd )

  call ncdf_fms_write ( Mgrid, Time, diag_id(pp_), gridm, ph )
  call ncdf_fms_write ( Mgrid, Time, diag_id(gz_), gridm, gz )
  call ncdf_fms_write ( Mgrid, Time, id_fz, gridm, fz )
  call ncdf_fms_write ( Mgrid, Time, diag_id(rh_), gridm, rh )
  call ncdf_fms_write ( Mgrid, Time, diag_id(sd_), gridm, sd )
  call ncdf_fms_write ( Mgrid, Time, diag_id(sm_), gridm, sm )
  call ncdf_fms_write ( Mgrid, Time, diag_id(om_), gridm, om )

!-----------------------------------------------------------------------
! write perturbation fields f(x,y,z) 
!-----------------------------------------------------------------------

  do k=kbd,ked
     do j=jbd,jed
        pp(:,j,k) = ph(:,j,k) - ppref(j,k)
        up(:,j,k) = uu(:,j,k) - uuref(j,k)
        tp(:,j,k) = ta(:,j,k) - taref(j,k)
        qp(:,j,k) = qv(:,j,k) - qvref(j,k)
        gp(:,j,k) = gz(:,j,k) - gzref(j,k)
     enddo
  enddo
gp = max(-15000.,min(15000.,gp))
  call ncdf_fms_write ( Mgrid, Time, pert_id(pp_), gridm, pp )
  call ncdf_fms_write ( Mgrid, Time, pert_id(uu_), gridu, up )
  call ncdf_fms_write ( Mgrid, Time, pert_id(ta_), gridm, tp )
  call ncdf_fms_write ( Mgrid, Time, pert_id(qv_), gridm, qp )
  call ncdf_fms_write ( Mgrid, Time, pert_id(qc_), gridm, gp )  ! use qc for gz

  call mpp_clock_end ( zetac_write_history_clock )

  if ( nsteps_stdout == 0 ) return
  n = n+1
!  if ( n > 1 .and. mod(n,nsteps_stdout) /= 0) return
  if (mod(n,nsteps_stdout) /= 0) return

!-----------------------------------------------------------------------
! surface diagnostics to stdout
!-----------------------------------------------------------------------

  call mpp_clock_begin ( zetac_history_stdout_clock )

  maxw = mxval*delp/100.
  write (stdout(), '(/,"sigma_history: time =",f10.0,2x,"(",f6.1,"% )",&
&                 6x,"maxval =",f6.2)')    elapsed_time, percent, maxw

  if (verbose) then

     call minmaxavg (                                                  &
                                                                Mgrid, &
                                                          ta(:,:,ked), &
                                                  hmask, area, varbar, &
                                       minlon, minlat, maxlon, maxlat, &
                                                           var1, var2 )

     write (stdout(),                                                  &
                '("sfc temp min/avg/max:",3f7.1," ... max at",2f6.1)') & 
                                     var1, varbar, var2, maxlon, maxlat

     call minmaxavg (                                                  &
                                                                Mgrid, &
                                                                   ps, &
                                                  hmask, area, varbar, &
                                       minlon, minlat, maxlon, maxlat, &
                                                           var1, var2 )

     write (stdout(),                                                  &
                '("sfc pres min/avg/max:",3f7.1," ... min at",2f6.1)') &
                     0.01*var1, 0.01*varbar, 0.01*var2, minlon, minlat

! energy tendencies

     call minmaxavg (                                                  &
                                                                Mgrid, &
                                                                  kdp, &
                                                    hmask, area, var1 )


     call minmaxavg (                                                  &
                                                                Mgrid, &
                                                                  pdp, &
                                                    hmask, area, var2 )

     varbar = var1 + var2
     detot = (varbar - esave)/(nsteps_stdout*0.5*delt)
     esave = varbar

     write (stdout(), '("total energy (10^6 J/m^2):",3f10.2)')         &
                                 1.e-6*var1, 1.e-6*var2, 1.e-6*varbar

     call minmaxavg (                                                  &
                                                                Mgrid, &
                                                                 ktdp, &
                                                    hmask, area, var1 )

     call minmaxavg (                                                  &
                                                                Mgrid, &
                                                                 ptdp, &
                                                    hmask, area, var2 )

     varbar = var1 + var2
     nsum = nsum+1
     if (nsum == 1) detot = varbar
     varsum = varsum + (varbar - detot)
     write (stdout(), '("energy tendencies (W/m^2):"4f10.1,f10.2)')    &
                               var1, var2, varbar, detot, varsum/nsum

! q tendency

!!$     call minmaxavg (                                                  &
!!$                                                                Mgrid, &
!!$                                                                 qtdp, &
!!$                                                    hmask, area, var2 )
!!$
!!$     dqtot = (qbar - qsave)/(nsteps_stdout*0.5*delt)
!!$     qsave = qbar
!!$
!!$     if (nsum ==  1) dqtot = var2
!!$     dqsum = dqsum + (var2 - dqtot)
!!$     write (stdout(),                                                  &
!!$             '("qvapor tendencies (10^-6 kg/m2/s):"12x,2f10.0,f10.2)') &
!!$                               1.e6*var2, 1.e6*dqtot, 1.e6*dqsum/nsum

     do j=jbc,jec
        do i=ibc,iec
           pdp(i,j) = minval(ps(i-1:i+1,j-1:j+1))
           p_is_min(i,j) = (ps(i,j) == pdp(i,j) .and. ps(i,j) < pcrit)
        enddo
     enddo
     var = 0.0
     var (mpp_pe()) = count(p_is_min(ibc:iec,jbc:jec))
     call mpp_sum ( var, size(var) )
     write (stdout(), '("vortices: ",i4)') int(sum(var))

  endif

  call mpp_clock_end ( zetac_history_stdout_clock )

  return
end subroutine write_history

!#######################################################################

subroutine write_static ( Mgrid, Time, ppref, uuref, taref, qvref )

type (horiz_grid_type),                         intent (in)    ::      &
                                                                Mgrid

type (time_type),                               intent (in)    ::      &
                                                                 Time

real, dimension(Mgrid%jbg:Mgrid%jeg, Mgrid%iz), intent (in)    ::      &
                                           ppref, uuref, taref, qvref

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

logical :: used

!-----------------------------------------------------------------------
! write static fields f(z)
!-----------------------------------------------------------------------

  if (grid_id(zetam_) > 0) used = send_data ( grid_id(zetam_),         & 
                                                          zetam, Time )
  if (grid_id(zetaw_) > 0) used = send_data ( grid_id(zetaw_),         & 
                                                          zetaw, Time )

!-----------------------------------------------------------------------
! write static fields f(y,z) - global axes (dangerous)
!-----------------------------------------------------------------------

  if (refs_id(ppref_) > 0) used = send_data ( refs_id(ppref_),         &
                                               ppref(1:iy,1:iz), Time )
  if (refs_id(uuref_) > 0) used = send_data ( refs_id(uuref_),         &
                                               uuref(1:iy,1:iz), Time )
  if (refs_id(taref_) > 0) used = send_data ( refs_id(taref_),         &
                                               taref(1:iy,1:iz), Time )
  if (refs_id(qvref_) > 0) used = send_data ( refs_id(qvref_),         &
                                               qvref(1:iy,1:iz), Time )

  return
end subroutine write_static

!#######################################################################

subroutine minmaxavg (                                                 &
                                                                Mgrid, &
                                                                field, &
                                                     hmax, atot, fbar, &
                                       minlon, minlat, maxlon, maxlat, &
                                                           fmin, fmax )

type (horiz_grid_type),                     intent (in)   ::           &
                                                                Mgrid

real, dimension(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed),              &
                                            intent (in)   ::           &
                                                                field

real, intent(in)  :: hmax, atot
real, intent(out) :: fbar
real, optional, intent(out) :: minlon, minlat, maxlon, maxlat
real, optional, intent(out) :: fmin, fmax

integer :: j, pe(1)
real, dimension(Mgrid%ibc:Mgrid%iec,Mgrid%jbc:Mgrid%jec) :: data

  data = field(ibc:iec,jbc:jec)

  do j=jbc,jec 
     tmp(j) = sum(data(:,j))!, mask=(topog(ibc:iec,j) <= hmax))*cosulat(j)
  enddo
  var = 0.0
  var (mpp_pe()) = sum(tmp)

  call mpp_sum ( var, size(var) )

  fbar = sum(var)/atot

  if ( present(fmin) ) then
     var = 0.0
     loc = 0
     var(mpp_pe())   = minval ( data, mask = (topog(ibc:iec,jbc:jec) <= hmax) )
     loc(:,mpp_pe()) = minloc ( data, mask = (topog(ibc:iec,jbc:jec) <= hmax) ) + offset
     call mpp_sum ( var, size(var) )
     call mpp_sum ( loc, size(loc) )
     pe = minloc(var)-1
     if ( present(minlon) ) minlon = vlon(loc(1,pe(1)))
     if ( present(minlat) ) minlat = ulat(loc(2,pe(1)))
     fmin = var(pe(1))
  endif

  if ( present(fmax) ) then
     var = 0.0
     loc = 0
     var(mpp_pe())   = maxval ( data, mask = (topog(ibc:iec,jbc:jec) <= hmax) )
     loc(:,mpp_pe()) = maxloc ( data, mask = (topog(ibc:iec,jbc:jec) <= hmax) ) + offset
     call mpp_sum ( var, size(var) )
     call mpp_sum ( loc, size(loc) )
     pe = maxloc(var)-1
     if ( present(maxlon) ) maxlon = vlon(loc(1,pe(1)))
     if ( present(maxlat) ) maxlat = ulat(loc(2,pe(1)))
     fmax = var(pe(1))
  endif

  return
end subroutine minmaxavg

!#######################################################################

subroutine write_history_init (                                        &
                                              Mgrid, Hmetric, Vmetric, &
                                                      Time, Time_step )

type (horiz_grid_type),                     intent (in)   ::           &
                                                                Mgrid

type (horiz_metric_type),                   intent (in)   ::           &
                                                              Hmetric

type (vert_metric_type),                    intent (in)   ::           &
                                                              Vmetric

type (time_type),                           intent (in)   ::           &
                                                      Time, Time_step

real    :: tolerance, tmin, tmax
integer :: unit, io, ierr, index
integer :: i, j, k
real :: ht

real, dimension(Mgrid%jbd:Mgrid%jed,Mgrid%iz) :: ppref

namelist /history_nml/ pcrit, nsteps_stdout, verbose, hmask

!-----------------------------------------------------------------------
! read namelist
!-----------------------------------------------------------------------

  read (input_nml_file, history_nml, iostat=io)
  ierr = check_nml_error(io,'history_nml')

   zetac_write_history_clock   = mpp_clock_id( 'zetac_write_history' )
   zetac_history_stdout_clock   = mpp_clock_id( 'zetac_history_stdout' )

!-----------------------------------------------------------------------
! write version number to logfile
!-----------------------------------------------------------------------

  call write_version_number (version, tag)

  ibc = Mgrid%ibc
  iec = Mgrid%iec
  ibd = Mgrid%ibd
  ied = Mgrid%ied

  jbc = Mgrid%jbc
  jec = Mgrid%jec
  jbd = Mgrid%jbd
  jed = Mgrid%jed

  kbd = Mgrid%kbd
  ked = Mgrid%ked

  iy = Mgrid%iy
  iz = Mgrid%iz

!-----------------------------------------------------------------------
! get height fields
!-----------------------------------------------------------------------

  allocate (zetaw(kbd:ked))
  allocate (zetam(kbd:ked))
  allocate (topog(ibd:ied,jbd:jed))

  zetam = Vmetric%zetam
  zetaw = Vmetric%zetaw
  topog = Vmetric%topog
  ptop  = Vmetric%ptop
  pbot  = Vmetric%pbot
  delp  = pbot - ptop

!-----------------------------------------------------------------------
! initialize axes
!-----------------------------------------------------------------------

  call axes_history ( Mgrid, Hmetric, Vmetric )

  axid_mm  = (/axid(xm_),axid(ym_)/)
  axid_uum = (/axid(xu_),axid(yu_),axid(zm_)/)
  axid_vvm = (/axid(xv_),axid(yv_),axid(zm_)/)
  axid_mmm = (/axid(xm_),axid(ym_),axid(zm_)/)
  axid_mmw = (/axid(xm_),axid(ym_),axid(zw_)/)
  axid_gm =  (/axid(yg_),axid(zm_)/)

!-----------------------------------------------------------------------
! register grid fields f(z)
!-----------------------------------------------------------------------

  call get_info ( grid_names(zetaw_), longname, units )
  grid_id(zetaw_) = register_static_field (                            &
              model, grid_names(zetaw_), axid_mmw(3:3),                &
              longname, units, silly, range=(/0.0,2.0/) )

  call get_info ( grid_names(zetam_), longname, units )
  grid_id(zetam_) = register_static_field (                            &
              model, grid_names(zetam_), axid_mmm(3:3),                &
              longname, units, silly, range=(/-1.0,2.0/) )

!-----------------------------------------------------------------------
! register static fields f(y,z)
!-----------------------------------------------------------------------

  call get_info ( refs_names(ppref_), longname, units )
  refs_id(ppref_) = register_static_field (                            &
              model, refs_names(ppref_), axid_gm,                      &
              longname, units, silly, range=(/-60.0,120.0/) )

  call get_info ( refs_names(uuref_), longname, units )
  refs_id(uuref_) = register_static_field (                            &
              model, refs_names(uuref_), axid_gm,                      &
              longname, units, silly, range=(/-60.0,120.0/) )

  call get_info ( refs_names(taref_), longname, units )
  refs_id(taref_) = register_static_field (                            &
              model, refs_names(taref_), axid_gm,                      &
              longname, units, silly, range=(/200.0,750.0/) )
 
  call get_info ( refs_names(qvref_), longname, units )
  refs_id(qvref_) = register_static_field (                            &
              model, refs_names(qvref_), axid_gm,                      &
              longname, units, silly, range=(/-0.005,0.03/) )

!-----------------------------------------------------------------------
! register prognostic fields f(x,y,z,t)
!-----------------------------------------------------------------------
 
  call get_info ( prog_names(uu_), longname, units )
  prog_id(uu_) = register_diag_field (                                 &
            model, prog_names(uu_), axid_uum, Time,                    &
            longname, units, silly, range=(/-120.0,120.0/) )

  call get_info ( prog_names(vv_), longname, units )
  prog_id(vv_) = register_diag_field (                                 &
            model, prog_names(vv_), axid_vvm, Time,                    &
            longname, units, silly, range=(/-120.0,120.0/) )
    
  call get_info ( prog_names(oo_), longname, units )
  prog_id(oo_) = register_diag_field (                                 &
            model, prog_names(oo_), axid_mmw, Time,                    &
            longname, units, silly, range=(/-40.0,40.0/) )
    
  call get_info ( prog_names(ta_), longname, units )
  prog_id(ta_) = register_diag_field (                                 &
            model, prog_names(ta_), axid_mmm, Time,                    &
            longname, units, silly, range=(/170.0,350.0/) )
    
  call get_info ( prog_names(qv_), longname, units )
  prog_id(qv_) = register_diag_field (                                 &
            model, prog_names(qv_), axid_mmm, Time,                    &
            longname, units, silly, range=(/-0.005,0.04/) )

  call get_info ( prog_names(qc_), longname, units )
  prog_id(qc_) = register_diag_field (                                 &
            model, prog_names(qc_), axid_mmm, Time,                    &
            longname, units, silly, range=(/-1.e5,1.e5/) )

!-----------------------------------------------------------------------
! register diagnostic fields f(x,y,z)
!-----------------------------------------------------------------------

  call get_info ( diag_names(pp_), longname, units )
  diag_id(pp_) = register_diag_field (                                 &
            model, diag_names(pp_), axid_mmw, Time,                    &
            longname, units, silly, range=(/0.0,110000.0/) )

  call get_info ( diag_names(gz_), longname, units )
  diag_id(gz_) = register_diag_field (                                 &
            model, diag_names(gz_), axid_mmw, Time,                    &
            longname, units, silly, range=(/0.0,2.6e5/) )

  name = "fz"
  longname = "geopotential at full-levels"
  units = "m2/s2"
  id_fz = register_diag_field (                                        &
            model, name, axid_mmm, Time,                               &
            longname, units, silly, range=(/0.0,4.0e5/) )

  call get_info ( diag_names(rh_), longname, units )
  diag_id(rh_) = register_diag_field (                                 &
            model, diag_names(rh_), axid_mmm, Time,                    &
            longname, units, silly, range=(/-0.05,2.00/) )

  call get_info ( diag_names(sd_), longname, units )
  diag_id(sd_) = register_diag_field (                                 &
            model, diag_names(sd_), axid_mmm, Time,                    &
            longname, units, silly, range=(/1.e5,1.e5/) )

  call get_info ( diag_names(sm_), longname, units )
  diag_id(sm_) = register_diag_field (                                 &
            model, diag_names(sm_), axid_mmm, Time,                    &
            longname, units, silly, range=(/1.e5,1.e5/) )

  call get_info ( diag_names(om_), longname, units )
  diag_id(om_) = register_diag_field (                                 &
            model, diag_names(om_), axid_mmm, Time,                    &
            longname, units, silly, range=(/-100.0,100.0/) )
    
  call get_info ( prog_names(uu_), longname, units )
  longname = trim(longname) // "at mass points"
  id_um = register_diag_field (                                        &
            model, 'um', axid_mmm, Time,                               &
            longname, units, silly, range=(/-130.0,130.0/) )

  call get_info ( prog_names(vv_), longname, units )
  longname = trim(longname) // "at mass points"
  id_vm = register_diag_field (                                        &
            model, 'vm', axid_mmm, Time,                               &
            longname, units, silly, range=(/-130.0,130.0/) )

  call get_info ( prog_names(uu_), longname, units )
  longname = "total velocity"
  id_vtot = register_diag_field (                                      &
            model, 'vtot', axid_mmm, Time,                             &
            longname, units, silly, range=(/-140.0,140.0/) )

  longname = "horizontal mass divergence"
  units = "kg s^-1 m^-2"
  id_hdiv = register_diag_field (                                      &
            model, 'hdiv', axid_mmm, Time,                             &
            longname, units, silly, range=(/-12.,12./) )

  longname = "vorticity"
  units = "s^-1"
  id_vort = register_diag_field (                                      &
            model, 'vort', axid_mmm, Time,                             &
            longname, units, silly, range=(/-0.1,0.1/) )

  longname = "kinetic energy per unit mass"
  units = 'm2/s2'
  id_ekin = register_diag_field (                                      &
            model, 'ekin', axid_mmm, Time,                             &
            longname, units, silly, range=(/-1.e9, 1.e9/) )

  longname = "potential energy per unit mass"
  units = 'm2/s2'
  id_epot = register_diag_field (                                      &
            model, 'epot', axid_mmm, Time,                             &
            longname, units, silly, range=(/-1.e9, 1.e9/) )

  longname = "column integrated kinetic energy"
  units = 'J/m2'
  id_ekin_int = register_diag_field (                                  &
            model, 'ekin_int', axid_mm, Time,                          &
            longname, units, silly, range=(/-1.e9, 1.e9/) )

  longname = "column integrated potential energy"
  units = 'J/m2'
  id_epot_int = register_diag_field (                                  &
            model, 'epot_int', axid_mm, Time,                          &
            longname, units, silly, range=(/-1.e9, 1.e9/) )

!-----------------------------------------------------------------------
! register perturbation fields f(x,y,z,t)
!-----------------------------------------------------------------------

  call get_info ( diag_names(pp_), longname, units )
  longname = "perturbation " // trim(longname)
  pert_id(pp_) = register_diag_field (                                 &
            model, trim(diag_names(pp_))//'_pert', axid_uum, Time,     &
            longname, units, silly, range=(/-75.0,75.0/) )
    
  call get_info ( prog_names(uu_), longname, units )
  longname = "perturbation " // trim(longname)
  pert_id(uu_) = register_diag_field (                                 &
            model, trim(prog_names(uu_))//'_pert', axid_uum, Time,     &
            longname, units, silly, range=(/-75.0,75.0/) )
    
  call get_info ( prog_names(ta_), longname, units )
  longname = "perturbation " // trim(longname)
  pert_id(ta_) = register_diag_field (                                 &
            model, trim(prog_names(ta_))//'_pert', axid_mmm, Time,     &
            longname, units, silly, range=(/-100.0,100.0/) )
    
  call get_info ( prog_names(qv_), longname, units )
  longname = "perturbation " // trim(longname)
  pert_id(qv_) = register_diag_field (                                 &
            model, trim(prog_names(qv_))//'_pert', axid_mmm, Time,     &
            longname, units, silly, range=(/-0.02,0.04/) )

  call get_info ( diag_names(gz_), longname, units )
  longname = "perturbation " // trim(longname)
  pert_id(qc_) = register_diag_field (                                 &
            model, trim(diag_names(gz_))//'_pert', axid_mmw, Time,     &
            longname, units, silly, range=(/-5.e4,5.e4/) )

!-----------------------------------------------------------------------
! energy fields
!-----------------------------------------------------------------------

  name = "kedt"
  longname = "kinetic energy generation"
  units = 'm2/s3'
  id_kedt = register_diag_field (                                      &
                      model, name, axid_mmm, Time,                     &
                      longname, units, silly, range=(/-20., 20./) )

  name = "pedt"
  longname = "potential energy generation"
  units = 'm2/s3'
  id_pedt = register_diag_field (                                      &
                      model, name, axid_mmm, Time,                     &
                      longname, units, silly, range=(/-50., 50./) )

  name = "kedt_int"
  longname = "column kinetic energy generation"
  units = 'W/m2'
  id_kedt_int = register_diag_field (                                  &
                      model, name, axid_mm, Time,                      &
                      longname, units, silly, range=(/-1.e6, 1.e6/) )

  name = "pedt_int"
  longname = "column potential energy generation"
  units = 'W/m2'
  id_pedt_int = register_diag_field (                                  &
                      model, name, axid_mm, Time,                      &
                      longname, units, silly, range=(/-1.e6, 1.e6/) )

!-----------------------------------------------------------------------
! allocate local arrays
!-----------------------------------------------------------------------

  allocate (udp(ibd:ied,jbd:jed))
  allocate (vdp(ibd:ied,jbd:jed))
  allocate (kdp(ibd:ied,jbd:jed))
  allocate (pdp(ibd:ied,jbd:jed))
  allocate (ktdp(ibd:ied,jbd:jed))
  allocate (ptdp(ibd:ied,jbd:jed))
  allocate (pp(ibd:ied,jbd:jed,kbd:ked))
  allocate (up(ibd:ied,jbd:jed,kbd:ked))
  allocate (tp(ibd:ied,jbd:jed,kbd:ked))
  allocate (qp(ibd:ied,jbd:jed,kbd:ked))
  allocate (gp(ibd:ied,jbd:jed,kbd:ked))
  allocate (ph(ibd:ied,jbd:jed,kbd:ked))
  allocate (pf(ibd:ied,jbd:jed,kbd:ked))
  allocate (sd(ibd:ied,jbd:jed,kbd:ked))
  allocate (vtot(ibd:ied,jbd:jed,kbd:ked))
  allocate (hdiv(ibd:ied,jbd:jed,kbd:ked))
  allocate (vort(ibd:ied,jbd:jed,kbd:ked))
  allocate (ekin(ibd:ied,jbd:jed,kbd:ked))
  allocate (epot(ibd:ied,jbd:jed,kbd:ked))
  allocate (ukin(ibd:ied,jbd:jed,kbd:ked))
  allocate (vkin(ibd:ied,jbd:jed,kbd:ked))
  allocate (epref(jbd:jed,kbd:ked))

  allocate (p_is_min(ibd:ied,jbc:jed))

  allocate (temp(0:mpp_npes()-1))
  allocate (temp1(0:mpp_npes()-1),temp2(0:mpp_npes()-1))
  allocate (var(0:mpp_npes()-1))
  allocate (loc(2,0:mpp_npes()-1))
  allocate (offset(2))

  offset = (/ibc-1, jbc-1/)

  allocate (vlon   (Mgrid%ibg:Mgrid%ieg))
  allocate (ulat   (Mgrid%jbg:Mgrid%jeg))
  allocate (cosulat(Mgrid%jbg:Mgrid%jeg))
  allocate (cosvlat(Mgrid%jbg:Mgrid%jeg))
  allocate (dxu    (Mgrid%jbg:Mgrid%jeg))

  ulat    = Hmetric%ulat
  vlon    = Hmetric%vlon
  cosulat = Hmetric%cosulat
  cosvlat = Hmetric%cosvlat
  dxu     = Hmetric%dxu
  dy      = Hmetric%dy

  allocate (tmp(jbc:jec))

  do j=jbc,jec
     tmp(j) = iec+1-ibc!count(topog(ibc:iec,j) <= hmask)*cosulat(j)
  enddo
  var = 0.0
  var (mpp_pe()) = sum(tmp)
  call mpp_sum ( var, size(var) )
  area = max( 0.1, sum(var) )
  ainf = Mgrid%ix*sum(cosulat(1:Mgrid%iy))

  if ( area == 0.0 ) then
     area = tiny
     if (mpp_pe() == 0) call error_mesg ( module,                      &
                            'averaging region has zero area', WARNING )
  endif

  call get_time ( Time_step, nsecs )

  return
end subroutine write_history_init

!#######################################################################

end module zetac_history_mod
