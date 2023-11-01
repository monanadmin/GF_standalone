# script to execute ${executable} sometimes with HPCToolkit profiler

# MONAN
#
# Author: Denis Eiras
#
# E-mail: denis.eiras@inpe.br  
#
# Date: 2023/02/08 
#
# Version: TODO -> USAR VERSÃO DA RELEASE, QUANDO FOR GERADA
#

# Full description
# =================
# This script executes locally ${executable} x times (configurable) with HPCToolkit or gprof profile. 
#
# Observe:
# - The execution must be done in the original folder (scripts)
#
# Usage: ./run_gf_times_with_hpctoolkit.sh <run_tool>
#
# Where: <run_tool> is 'hpctoolkit' or 'gprof'
#
#
# History
# =======
# 2023/02/08 - Initial commit - deniseiras
# 2023/02/16 - Hpctoolkit and gprof code
# 2023/11/01 - Including parameter options and gprof submission job
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
run_tool=$1
times_to_execute=1

if [ -z ${run_tool} ]; then
  echo "No argument providade. Please, inform 'hpctoolkit' or 'gprof'"
  exit -1
elif [ "$1" != "hpctoolkit" ] && [ "$1" != "gprof" ]; then
    echo "Please, inform 'hpctoolkit' or 'gprof' in the first argument. Informed: $1"
fi


source env_hpctoolkit_gprof.sh

# code
ln -fs ../datain/gf.inp
ln -fs ../datain/GF_ConvPar_nml
# ln -fs ../datain/GATE.dat
ln -fs ../datain/gf_dataIn-9600.bin

executable="gf.x"
exec_patch="../bin/gf.x"

hpct_measure="hpctoolkit-${executable}-measurements"
hpct_struct="${executable}.hpcstruct"
hpct_db="hpctoolkit-${executable}-database"

NOW=$(date +'%y%m%d-%H%M%S')
mv executions executions_$NOW

for counter in $(seq 1 $times_to_execute); do 
  dir_exec="executions/exec_${counter}"
  rm -rf ${dir_exec}
  mkdir -p ${dir_exec}

  echo -e "\n\n\nExecution of profiler ${run_tool} on ${executable} # $counter"

  if [ ${run_tool} == 'hpctoolkit' ]; then
    # for HPCtoolkit on nodes  =======================
    sbatch -W submit_hpctoolkit_from_script_run_gf_times.sbatch
    sleep 30
    mv ${hpct_measure} ${hpct_struct} ${hpct_db} ${dir_exec}
    echo -e "\n\n\nFinished!!! \n"
    echo -e "\n\nCheck HPCToolkit results:\n\tsource env_hpctoolkit_gprof.sh\n\thpcviewer ${dir_exec}/hpctoolkit-${executable}-database"

  elif [ ${run_tool} == 'gprof' ]; then
    # for gprof on nodes  =======================
    sbatch -W submit_gprof_from_script_run_gf_times.sbatch
    sleep 30
    gprof --graph ${exec_path} > gprof.out  # other options: --exec-times --graph --brief --flat-profile
    sleep 5
    mv ./gmon.out ./gprof.out ${dir_exec}
    echo -e "\n\n\nFinished!!! \n"
    echo -e "\n\nCheck Gprof results     :\n\tmore ${dir_exec}/gprof.out"
  fi
  
  #mv ref_g.* ${dir_exec}
  mv slurm*.out ${dir_exec}
  echo -e "\n\n"

done



# for running locally

  # for HPCtoolkit local  =======================
  # hpcrun -t ../bin/${executable}
  # hpcstruct ../bin/${executable}
  # hpcprof -I . -S ${hpct_struct} ${hpct_measure}
  # mv $hpct_measure $hpct_structf ${hpct_db} ${dir_exec}

  # for gprof local =============================
  # ../bin/${executable}
  # gprof --graph ../bin/${executable} > gprof.out  # --exec-times --graph --brief --flat-profile
  # mv ./gmon.out ./gprof.out ${dir_exec}

