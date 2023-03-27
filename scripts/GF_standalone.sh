#!/bin/sh

# Script que compila e executa o codigo GF_StandAlone:

#if [ $# -ne 1 ]
#then
#   echo ""
#   echo "/scripts/clima/tools/makeproj [PROJECT root name]"
#   echo ""
#   exit
#fi
 
DIRHOME=$PWD/..
SCRIPTS=${DIRHOME}/scripts
DATAOUT=${DIRHOME}/dataout
DATAIN=${DIRHOME}/datain
SRC=${DIRHOME}/src
BIN=${DIRHOME}/bin


# Verificando o argumento de entrada
if [ -z "${1}" ]
then
  echo "COMPILER is not set: gnu or intel" 
  exit 1
else
  COMPILER=`echo $1 | awk '{print $1}' `
  echo "COMPILER=$COMPILER"
fi


cd ${BIN}
rm -f ${BIN}/gf.x
echo "Compilando"
make clean
make $COMPILER
sleep 3


cd ${DATAIN}
echo "Executando"
time ${BIN}/gf.x

echo "done"
