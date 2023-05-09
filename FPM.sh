# bidu@lncc.br, abril 2023, proj monan
# uso do arquivo projectOpt.rsp 
# para escolha de parametros para o fpm
#   foi a unica forma encontrada de escolher os flags da compilação

rm build -r
#echo y | fpm  clean  
fpm  build  @projectOpt
comando="find build -name gf_teste |xargs -I{} cp \"{}\" bin/gf_teste"
echo $comando; eval $comando
