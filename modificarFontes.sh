#!/bin/bash
#plano para  fazer a alteração codigo para excluir os 
#   if + cycle

#134 locais encontrados para fazer as substituições 
#
#onde houver   : do i = its, itf
#substituir por: do vtp_index = 1, get_num_elements()
#
#
#onde houver : if (ierr(i) /= 0) cycle
#substituir por: i=get_index_value(vtp_index)
#
#plano:
#1. encontrar os numeros das linhas que possuem:
#    do i = its, itf
#   e colocar em um arquivo linhasOcorr.txt
#
fOrig=src/00modConvParGF.B90
fModif=src/modConvParGF.F90
cp $fOrig $fModif
padrao="do"; sep="-"
grep -n -E "ierr.*cycle"  src/modConvParGF.F90 -B 9|\
        grep -E "its, *itf|ierr.*\) cycle" |\
        grep -v -E "applyEnvSubs|Nunca| *! |or|1745" |\
        grep $padrao |awk -F${sep} '{print $1}' > linhasOcorrDO.txt

for nl in $(cat linhasOcorrDO.txt |head -99 )
  do
          novaLinha=" do vtp_index = 1, get_num_elements() ; i=get_index_value(vtp_index) !BD_n"
          comandoSed="${nl}s/do i.*itf/${novaLinha}/"
          comando="sed -i    '$comandoSed'  src/modConvParGF.F90"
	  echo $comando; 
	  eval $comando;
  done

padrao="if"; sep=":"
grep -n -E "ierr.*cycle"  src/modConvParGF.F90 -B 9|\
        grep -E "its, *itf|ierr.*\) cycle" |\
        grep -v -E "applyEnvSubs|Nunca| *! |or|1745" |\
        grep $padrao |awk -F${sep} '{print $1}' > linhasOcorrIF.txt

for nl in $(cat linhasOcorrIF.txt |head -99 )
  do
          comandoSed="${nl}s/^/!BD_n/"
          comando="sed -i    '$comandoSed'  src/modConvParGF.F90"
	  echo $comando; 
	  eval $comando;
done

