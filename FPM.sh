# bidu@lncc.br, abril 2023, proj monan
# para escolha de parametros para o fpm
# uso do arquivo projectOpt.rsp 
# ou por parametro na linha de comando

EXE="gf_teste"
EXEF="gf_teste.x"
echo Criando o executavel 
comando="echo y | fpm  clean; rm build/cache.*  "
echo $comando; eval $comando
comando="fpm  build  @projectOpt"
comando="fpm  build  -flag=\"-fPIC -O0 -g -fcheck=bounds -fbacktrace -ffree-line-length-200\""
echo $comando; eval $comando
comando="find build -name $EXE |xargs -I{} cp \"{}\" bin/$EXEF"
EXE="$(find build -name $EXE|head -n 1)"
#echo $EXE; exit
echo copiando $EXE para bin/$EXEF 
comando="cp $EXE bin/$EXEF"
echo $comando; eval $comando
comando="ls -ltr bin/$EXEF"
echo $comando; eval $comando
