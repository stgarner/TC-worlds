program zetac_coldstart

use fms_mod,    only : fms_init, fms_end,                              &
                       write_version_number, stdout
use mpp_mod,    only : mpp_init, mpp_pe
use fms_io_mod, only : fms_io_exit

use diag_manager_mod, only : diag_manager_init
use time_manager_mod, only : time_type, time_manager_init
use constants_mod,    only : constants_init

use zetac_axes_mod,              only : gridu, gridv, gridm
use zetac_domains_mod,           only : get_domains
use zetac_horiz_grid_type_mod,   only : horiz_grid_type
use zetac_horiz_metric_mod,      only : horiz_metric
use zetac_horiz_metric_type_mod, only : horiz_metric_type
use zetac_namelist_mod,          only : namelist_read
use zetac_set_topog_mod,         only : set_topog
use zetac_times_mod,             only : get_times
use zetac_update_halos_mod,      only : update_halos, update_halos_init
use zetac_user_perturb_mod,      only : user_perturb
use zetac_user_ref_atmos_mod,    only : user_ref_atmos
use zetac_vert_metric_mod,       only : vert_metric
use zetac_vert_metric_type_mod,  only : vert_metric_type
use zetac_write_coldstart_mod,   only : write_coldstart,               &
                                        coldstart_init, coldstart_end

implicit none

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

type (horiz_grid_type)   :: Mgrid
type (horiz_metric_type) :: Hmetric
type (vert_metric_type)  :: Vmetric
type (time_type)         :: Time_init, Time_end

real, allocatable, dimension(:)     :: psref, tsref, zetaw

real, allocatable, dimension(:,:)   :: uuref, taref, rhref

real, allocatable, dimension(:,:)   :: ppert

real, allocatable, dimension(:,:,:) :: upert, vpert,tpert, qpert, hpert 

real, allocatable, dimension(:,:)   :: topog, sst

integer :: ibd, ied, ix
integer :: jbd, jed, iy
integer :: kbd, ked, iz
logical :: lxopen, lyopen
logical :: do_cartesian=.false., do_cylindrical=.false.
integer :: nsteps, nbuf=2
real    :: rlonmin, rlonmax, rlatmin, rlatmax
real    :: pbot, ptop
real    :: zeta(20)

!-----------------------------------------------------------------------
! version control information
!-----------------------------------------------------------------------

character(len=128)  :: version =  '$Id: zetac_coldstart.f90,v 1.1.2.6.2.9 2005/08/07 00:13:25 stg Exp $'
character(len=128)  :: tag     =  '$Name:  $'

!-----------------------------------------------------------------------
! initialize mpp and fms
!-----------------------------------------------------------------------
 
  call mpp_init()

  call fms_init

!-----------------------------------------------------------------------
! write version number to logfile
!-----------------------------------------------------------------------

  call write_version_number (version, tag)

!-----------------------------------------------------------------------
! initialize time variables
!-----------------------------------------------------------------------

  call get_times ( Time_init, Time_end )

!-----------------------------------------------------------------------
! read remining namelists
!-----------------------------------------------------------------------

  call namelist_read (                                                 &
                                                           ix, iy, iz, &
                                                         nsteps, nbuf, &
                         lxopen, lyopen, do_cartesian, do_cylindrical, &
                                   rlonmin, rlonmax, rlatmin, rlatmax, &
                                                     pbot, ptop, zeta )

  allocate (zetaw(iz))
  zetaw = zeta(1:iz)

!-----------------------------------------------------------------------
! get domain decomposition and define associated derived type
!-----------------------------------------------------------------------

  call get_domains (                                                   &
                                                     ix, iy, iz, nbuf, &
                                                Mgrid, lxopen, lyopen )

!-----------------------------------------------------------------------
! initialize wrapper for mpp_update_domains
!-----------------------------------------------------------------------

  call update_halos_init ( Mgrid )

!-----------------------------------------------------------------------
! iniitialize diagnostics manager (needed for base date)
!-----------------------------------------------------------------------

  call diag_manager_init

!-----------------------------------------------------------------------
! initialize constants
!-----------------------------------------------------------------------

  call constants_init

!-----------------------------------------------------------------------
! initialize time_manager
!-----------------------------------------------------------------------

  call time_manager_init

!-----------------------------------------------------------------------
! define horizontal metric terms and associated derived type variable
!-----------------------------------------------------------------------

  call horiz_metric (                                                  &
                                                       Mgrid, Hmetric, &
                                   rlonmin, rlonmax, rlatmin, rlatmax, &
                                         do_cartesian, do_cylindrical )

!-----------------------------------------------------------------------
! initialize topography and sst
!-----------------------------------------------------------------------

  allocate (topog(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed))
  allocate (  sst(Mgrid%ibd:Mgrid%ied,Mgrid%jbd:Mgrid%jed))

  call set_topog ( Mgrid, Hmetric, topog, sst )

!-----------------------------------------------------------------------
! define vertical metric terms and associated derived type variable
!-----------------------------------------------------------------------

  call vert_metric (                                                   &
                                              Mgrid, Hmetric, Vmetric, &
                                                           pbot, ptop, &
                                                          zetaw,topog )

!-----------------------------------------------------------------------
! initialize coldstart routine
!-----------------------------------------------------------------------

  call coldstart_init (                                                &
                                              Mgrid, Hmetric, Vmetric, &
                                   rlonmin, rlonmax, rlatmin, rlatmax, &
                                         do_cartesian, do_cylindrical )

!-----------------------------------------------------------------------
! allocate 2D reference atmosphere
!-----------------------------------------------------------------------

  allocate (psref(Mgrid%jbg:Mgrid%jeg))
  allocate (tsref(Mgrid%jbg:Mgrid%jeg))
  allocate (uuref(Mgrid%jbg:Mgrid%jeg,iz))
  allocate (taref(Mgrid%jbg:Mgrid%jeg,iz))
  allocate (rhref(Mgrid%jbg:Mgrid%jeg,iz))

!-----------------------------------------------------------------------
! define 2D reference atmosphere
!-----------------------------------------------------------------------

  call user_ref_atmos (                                                &
                                              Mgrid, Hmetric, Vmetric, &
                                    psref, tsref, uuref, taref, rhref )

!-----------------------------------------------------------------------
! allocate 3D reference atmosphere
!-----------------------------------------------------------------------

  ibd = Mgrid%ibd
  ied = Mgrid%ied
  jbd = Mgrid%jbd
  jed = Mgrid%jed
  kbd = Mgrid%kbd
  ked = Mgrid%ked

!-----------------------------------------------------------------------
! allocate 3D initial perturbation
!-----------------------------------------------------------------------

  allocate (ppert(ibd:ied,jbd:jed))
  allocate (upert(ibd:ied,jbd:jed,kbd:ked))
  allocate (vpert(ibd:ied,jbd:jed,kbd:ked))
  allocate (tpert(ibd:ied,jbd:jed,kbd:ked))
  allocate (hpert(ibd:ied,jbd:jed,kbd:ked))

!-----------------------------------------------------------------------
! define 3D perturbation
!-----------------------------------------------------------------------

  call user_perturb (                                                  &
                                              Mgrid, Hmetric, Vmetric, &
                                                  psref, taref, rhref, &
                                    ppert, upert, vpert, tpert, hpert )

!-----------------------------------------------------------------------
! interpolate fields to zeta surfaces
!-----------------------------------------------------------------------

  call update_halos ( Mgrid, ppert )
  call update_halos ( Mgrid, upert )
  call update_halos ( Mgrid, vpert )
  call update_halos ( Mgrid, tpert )
  call update_halos ( Mgrid, hpert )

!-----------------------------------------------------------------------
! process fields and write coldstart file
!-----------------------------------------------------------------------

  call write_coldstart (                                               &
                                            Mgrid, Vmetric, Time_init, &
                                                           topog, sst, &
                                    psref, tsref, uuref, taref, rhref, &
                                    ppert, upert, vpert, tpert, hpert )

!-----------------------------------------------------------------------
! exit coldstart routine
!-----------------------------------------------------------------------

  call coldstart_end

!-----------------------------------------------------------------------
! exit fms
!-----------------------------------------------------------------------

  call fms_io_exit

  call fms_end

  stop
end program zetac_coldstart
