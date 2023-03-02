#!/bin/sh

# Script que compila e executa o codigo GF_StandAlone:

#if [ $# -ne 1 ]
#then
#   echo ""
#   echo "/scripts/clima/tools/makeproj [PROJECT root name]"
#   echo ""
#   exit
#fi
 

DIRHOME=/mnt/beegfs/carlos.souza/GF_standalone
SCRIPTS=${DIRHOME}/scripts
DATAOUT=${DIRHOME}/dataout
DATAIN=${DIRHOME}/datain
SRC=${DIRHOME}/src
BIN=${DIRHOME}/bin

cd ${BIN}


echo "Compilando"
#make clean
rm -f gf.x
make -f Makefile_r8

echo "Executando"
time ./gf.x

echo "done"
