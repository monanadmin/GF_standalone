#!/bin/sh

# Script que compila e executa o codigo GF_StandAlone:

#if [ $# -ne 1 ]
#then
#   echo ""
#   echo "/scripts/clima/tools/makeproj [PROJECT root name]"
#   echo ""
#   exit
#fi
 

DIRHOME=$(pwd)/..
SCRIPTS=${DIRHOME}/scripts
DATAOUT=${DIRHOME}/dataout
DATAIN=${DIRHOME}/datain
SRC=${DIRHOME}/src
BIN=${DIRHOME}/bin


cd ${BIN}

echo "Compilando"
#make clean
make 
sleep 3


cd ${DATAIN}
echo "Executando"
time ${BIN}/gf.x

echo "done"
