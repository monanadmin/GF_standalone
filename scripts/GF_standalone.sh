#!/bin/sh

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
if [ -z "${1}" ]
then
  echo "Compiler is not set: gnu or intel"
  echo "gnu is set by default" 
  COMPILER='gnu'
else
  COMPILER=`echo $1 | awk '{print $1}' `
fi
  
echo "COMPILER=$COMPILER"

cd ${BIN}
#rm -f gf.x
rm -f gf_testex
echo "Compilando"
#make clean; make $COMPILER
(cd ../; ./FPM.sh)
sleep 3


cd ${DATAIN}
echo "Executando"
time ${BIN}/gf.x

echo "done"
