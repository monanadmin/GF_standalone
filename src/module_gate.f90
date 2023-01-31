module module_gate
   implicit none

   !-for BRAMS runs, set use_gate=.false.
   logical, parameter :: use_gate = .true.

   !- Here are the place for data related with the GATE soundings
   integer, parameter :: gate = 1   ! flag to turn on/off : 1/0
   integer, parameter :: klon = 161 ! number of soundings for gate
   integer, parameter :: klev = 41  ! number of vertical levels
   integer, parameter :: ktrac = 2   ! number of chemical tracers
   integer, parameter :: levs = klev
   integer, parameter :: nvar_grads = 300

   type cupout_vars
      real, pointer     :: varp(:, :)
      character(LEN=80), allocatable :: varn(:)
   end type cupout_vars

   type(cupout_vars), allocatable :: cupout(:)

   real, dimension(klon, klev):: pgeo, ppres, ptemp, pq, pu, pv, pvervel, &
                                 zrvten, ztten, zq1, zq2, zqr, zadvt, zadvq

   integer :: JL, KLEV_SOUND
   character(len=128) :: runname, runlabel, rundata = "NONE"

end module module_gate
