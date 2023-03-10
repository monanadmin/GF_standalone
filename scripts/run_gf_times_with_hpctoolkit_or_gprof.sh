# script to execute gf.x sometimes with HPCToolkit profiler

# MONAN
#
# Author: Denis Eiras
#
# E-mail: denis.eiras@inpe.br  
#
# Date: 2023/02/08 
#
# Version: TODO -> numero_da_versão_em_git_flow
#

# Full description
# =================
# This script executes locally gf.x 3 times (configurable) with HPCToolkit profile on. 
#
# Observe:
# - Must load HPCToolkit module before. 
#   ex: module load hpctoolkit-2021.10.15-gcc-9.4.0-i3utntd
# - The execution must be done in the original folder (scripts)
#
# Usage: ./run_gf_times_with_hpctoolkit.sh
#
#
# History
# =======
# 2023/02/08 - Initial commit - deniseiras
# 2023/02/16 - Inseridas opções para trabalhar com gprof e com Hpctoolkit no slurm
#
#
# Licence
# =======
#
# <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the  Free  Software  Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it  will be useful, but
# ** WITHOUT  ANY  WARRANTY **;  without  even  the   implied   warranty  of
# **MERCHANTABILITY** or **FITNESS FOR A  PARTICULAR PURPOSE**.  See  the, GNU
# GNU General Public License for more details.
#
# You should have received a copy  of the GNU General  Public  License
# along with this program.  If not, see [GNU Public    License](https://www.gnu.org/licenses/gpl-3.0.html).
#

# parameters
times_to_execute=1

source env_hpctoolkit.sh

# code
ln -fs ../datain/gf.inp
ln -fs ../datain/GF_ConvPar_nml
ln -fs ../datain/GATE.dat

hpct_measure="hpctoolkit-gf.x-measurements"
hpct_struct="gf.x.hpcstruct"
hpct_db="hpctoolkit-gf.x-database"

rm -rf executions 

for counter in $(seq 1 $times_to_execute); do 
  dir_exec="executions/exec_${counter}"
  rm -rf ${dir_exec}
  mkdir -p ${dir_exec}

  echo -e "\n\n\nExecution of profiler on gf.x # $counter"

  # for HPCtoolkit on nodes  =======================
  sbatch -W submit_hpctoolkit_from_script_run_gf_times.sbatch
  mv $hpct_measure $hpct_structf ${hpct_db} slurm*.out ${dir_exec}
    
  # for HPCtoolkit local  =======================
  # hpcrun -t ../bin/gf.x 
  # hpcstruct ../bin/gf.x 
  # hpcprof -I . -S ${hpct_struct} ${hpct_measure}
  # mv $hpct_measure $hpct_structf ${hpct_db} ${dir_exec}

  # for gprof local =============================
  # ../bin/gf.x
  # gprof --graph ../bin/gf.x > gprof.out  # --exec-times --graph --brief --flat-profile
  # mv ./gmon.out ./gprof.out ${dir_exec}

  # for all
  mv ref_g.* ${dir_exec}

done


echo -e "\n\n\nFinished!!! \nCheck results in executions folder, which could be databases of HPCToolkit and/or gprof.out, depending of code comented or not in the script"
echo -e "\n\nView HPCToolkit results:\n\tsource env_hpctoolkit.sh\n\thpcviewer ${dir_exec}/hpctoolkit-gf.x-database"
echo -e "\n\nView Gprof results     :\n\tmore ${dir_exec}/gprof.out"
echo -e "\n\n"
