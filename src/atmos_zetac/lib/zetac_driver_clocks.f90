module zetac_driver_clocks_mod

!-----------------------------------------------------------------------
! module information
!-----------------------------------------------------------------------

use mpp_mod, only : mpp_clock_id, MPP_CLOCK_DETAILED

implicit none

integer :: id_check_blowup
integer :: id_moisture_driver
integer :: id_time_filter
integer :: id_update_domains
integer :: id_history
integer :: id_write_restart
integer :: id_zetac_init

integer :: id_zetac_initialize
integer :: id_zetac_time_loop

integer :: flags = 0

!-----------------------------------------------------------------------
! version control information
!-----------------------------------------------------------------------

contains

!#######################################################################

subroutine driver_clocks_init

! flags = MPP_CLOCK_DETAILED

  id_check_blowup      = mpp_clock_id( 'check_blowup'      , flags=flags, grain = 20 )
  id_moisture_driver   = mpp_clock_id( 'moisture_driver'   , flags=flags, grain = 20 )
  id_time_filter       = mpp_clock_id( 'time_filter'       , flags=flags, grain = 20 )
  id_update_domains    = mpp_clock_id( 'update_domains'    , flags=flags, grain = 20 )
  id_history           = mpp_clock_id( 'history'           , flags=flags, grain = 20 )
  id_write_restart     = mpp_clock_id( 'write_restart',      flags=flags, grain = 20 )
  id_zetac_init        = mpp_clock_id( 'zetac_init'        , flags=flags, grain = 20 )

  id_zetac_initialize  = mpp_clock_id( 'initialization'    , flags=flags, grain = 10 )
  id_zetac_time_loop   = mpp_clock_id( 'time loop'         , flags=flags, grain = 10 )
  
end subroutine driver_clocks_init

!#######################################################################

end module zetac_driver_clocks_mod
