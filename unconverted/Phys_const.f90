module Phys_const

  implicit none

  real, parameter ::      &
       rgas    = 287.,    &
       cp      = 1004.,   &
       rm      = 461.,    &
       p00     = 1.e5,    &
       tcrit   = 273.15,  &
       g       = 9.80,    &
       cpor    = cp/rgas, &
       PKDCUT  = 75.,     &
       day_sec = 86400,   &
       xl      = 2.5e6

end module Phys_const
