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
if [ -z "${1}" ]
then
  echo "Compiler is not set: gnu or intel"
  echo "$COMPILER is set by default" 
fi
  
echo "COMPILER=$COMPILER"

cd ${BIN}
rm -f gf.x
echo "Compilando"
comando="make clean; make $COMPILER"
echo $comando; eval $comando

#(cd ../; ./FPM.sh)

#rm ${DATAOUT}/ref_g.gra
rm ${DATAOUT}/*.gra
cd ${DATAIN}
echo "Executando"
time ${BIN}/gf.x

echo "done"
