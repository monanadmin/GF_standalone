module modGate
   !! brief
   !!
   !! @note
   !!
   !! **Project**: MONAN
   !! **Author(s)**: Rodrigues, L.F. [LFR]
   !! **e-mail**: <mailto:luiz.rodrigues@inpe.br>
   !! **Date**:  03Fevereiro2023 11:03
   !!
   !! **Full description**:
   !! brief
   !!
   !! @endnote
   !!
   !! @warning
   !!
   !!  [](https://www.gnu.org/graphics/gplv3-127x51.png'')
   !!
   !!     This program is free software: you can redistribute it and/or modify
   !!     it under the terms of the GNU General Public License as published by
   !!     the  Free  Software  Foundation, either version 3 of the License, or
   !!     (at your option) any later version.
   !!
   !!     This program is distributed in the hope that it  will be useful, but
   !!     WITHOUT  ANY  WARRANTY;  without  even  the   implied   warranty  of
   !!     MERCHANTABILITY or FITNESS FOR A  PARTICULAR PURPOSE.  See  the, GNU
   !!     GNU General Public License for more details.
   !!
   !!     You should have received a copy  of the GNU General  Public  License
   !!     along with this program.  If not, see <https://www.gnu.org/licenses/>.
   !!
   !! @endwarning

   implicit none

   character(len=*), parameter :: sourceName = 'module_gate.f90' ! Nome do arquivo fonte
   character(len=*), parameter :: moduleName = 'modGate' ! Nome do módulo

   !-for BRAMS runs, set use_gate=.false.
   logical, parameter :: p_use_gate = .true.

   !- Here are the place for data related with the GATE soundings
   integer, parameter :: p_gate = 1   ! flag to turn on/off : 1/0
   integer, parameter :: p_klon = 161 ! number of soundings for gate
   integer, parameter :: p_klev = 41  ! number of vertical levels
   integer, parameter :: p_ktrac = 2   ! number of chemical tracers
   integer, parameter :: p_levs = p_klev
   integer, parameter :: p_nvar_grads = 300

   type t_cupout_vars
      real, pointer :: varp(:, :)
      character(len=80), allocatable :: varn(:)
   end type t_cupout_vars

   type(t_cupout_vars), allocatable :: cupout(:)

   real, dimension(p_klon, p_klev):: pgeo, ppres, ptemp, pq, pu, pv, pvervel, &
                                 zrvten, ztten, zq1, zq2, zqr, zadvt, zadvq

   integer :: jl, klev_sound
   character(len=128) :: runname, runlabel, rundata = "NONE"

contains

end module modGate
