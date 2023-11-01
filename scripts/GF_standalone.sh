#!/bin/bash

# Script que compila e executa o codigo GF_StandAlone:

# Usage: 
#
#        ./GF_standalone.sh gnu    :: para compilar com gfortran (default)
#        ./GF_standalone.sh intel  :: para compilar com ifort


DIRHOME=$PWD/..
SCRIPTS=${DIRHOME}/scripts
DATAOUT=${DIRHOME}/dataout
DATAIN=${DIRHOME}/datain
SRC=${DIRHOME}/src
BIN=${DIRHOME}/bin

# Verificando o argumento de entrada
COMPILER=${1:-"gnu"}
EXECUTE_LOCAL=${2:-true}

if [ -z "${1}" ]
then
  echo "Compiler is not set: gnu or intel"
  echo "$COMPILER is set by default" 
fi
echo "COMPILER=$COMPILER"

. env_hpctoolkit_gprof.sh

cd ${BIN}
/bin/cp Makefile_3D Makefile
echo "Compilando"
comando="make clean; make $COMPILER"
#comando="make $COMPILER"
echo $comando; eval $comando

if [ "${EXECUTE_LOCAL}" = true ]; then

  # for FPM Fortran Benchmarking use only
  #(cd ../; ./FPM.sh)
  
  (cd ${DATAOUT}; rm *.gra *.ctl)
  cd ${DATAIN}
  echo "Executando"
  time ${BIN}/gf.x

  echo "done"
fi
