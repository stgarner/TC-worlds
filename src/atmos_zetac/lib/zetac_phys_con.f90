module zetac_phys_con_mod

use constants_mod, only : kappa, rdgas, rvgas, cp_air

implicit none

real, parameter ::             &
   pref = 1013.e2,             &
   tref = 255.

real, parameter ::             &
   rrat   = rvgas/rdgas-1.0,   &  !  virtual effect on density
   eps    = rdgas/rvgas           !  R/Rv

real, parameter ::             &
   cp_H2O = 1880.00,           &  !  heat capacity of water vapor
   crat   = cp_H2O/cp_air-1.0

real, parameter ::             &
   tnuke  =   20.00               !  temperature offset for ice nucleation [K]

end module zetac_phys_con_mod
