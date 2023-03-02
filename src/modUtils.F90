module modUtils
   !! ## Modules with a lot of utilities
   !!
   !! ![](https://i.ibb.co/LNqGy3S/logo-Monan-Color-75x75.png)
   !! ## MONAN
   !!
   !! Author: Rodrigues, L.F. [LFR]
   !!
   !! E-mail: <mailto:luizfrodrigues@protonmail.com>
   !!
   !! Date: 02Março2023 16:55
   !!
   !! #####Version: 0.1.0
   !!
   !! ---
   !! **Full description**:
   !!
   !! Modules with a lot of utilities
   !!
   !! ** History**:
   !!
   !! --- 
   !! ** Licence **:
   !!
   !!  <img src="https://www.gnu.org/graphics/gplv3-127x51.png width="63">
   !!
   !!  This program is free software: you can redistribute it and/or modify
   !!  it under the terms of the GNU General Public License as published by
   !!  the  Free  Software  Foundation, either version 3 of the License, or
   !!  (at your option) any later version.
   !!
   !!  This program is distributed in the hope that it  will be useful, but
   !!  ** WITHOUT  ANY  WARRANTY **;  without  even  the   implied   warranty  of
   !!  **MERCHANTABILITY** or **FITNESS FOR A  PARTICULAR PURPOSE**.  See  the, GNU
   !!  GNU General Public License for more details.
   !!
   !!  You should have received a copy  of the GNU General  Public  License
   !!  along with this program.  If not, see [GNU Public License](https://www.gnu.org/licenses/gpl-3.0.html).
   !!
#ifdef MPI
   use mpi
#endif
   implicit none
   !use modConstants, only: only_list
   character(len=*), parameter :: p_source_name = 'modUtils.F90' 
   !! Source code name 
   character(len=*), parameter :: p_module_name = 'modUtils' 
   !! module name 

   private
   public :: StopExecution

contains

   ! --------------------------------------------------------------------- 
   function StopExecution(message, source, proced, processor) result(end_status)
        !! ## Para uma execução e imprime uma mensagem
        !!
        !! Author: Rodrigues, L.F. [LFR]
        !! 
        !! e-mail: <mailto:luiz.rodrigues@inpe.br>
        !!
        !! Date:  13Outubro2022 11:03
        !!
        !! —
        !! **Full description**:
        !!
        !! Stop a execution and print a message
        !!
        !! —
        implicit none
        ! Parameters:
        character(len=*), parameter :: p_procedure_name = 'StopExecution' 
     
        ! Variables (input):
        character(len=*), intent(in) :: message
        !! Mensagem a ser impressa
        character(len=*), intent(in) :: source
        !! Nome do arquivo fonte
        character(len=*), intent(in) :: proced
        !! Nome do procedure
        integer, intent(in), optional :: processor
        !! Se é paralelizada o número do processador deve estar presente
     
        ! Local variables:
        integer :: ierr
        integer :: end_status
     
        end_status = 0
        if(present(processor)) then
#ifdef MPI
           call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
           if(processor == 0) then
              print *,'##############################################################'
              print *, "Fortran Stop - "//source//" : "//proced//"()"
              print *, message
              print *,'##############################################################'
              print *,''
#ifdef MPI
              call mpi_finalize(ierr)
#else
              stop
#endif
           endif
        else
           print *,'####################################################################'
           print *, "Fortran Stop - "//source//" : "//proced//"()"
           print *, message
           print *,'####################################################################'
           print *,''
           stop 
        endif
   end function StopExecution


end module modUtils