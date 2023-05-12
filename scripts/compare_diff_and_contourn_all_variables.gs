
* Script to comparte grads variables between 2 files
* The script iterates over levels, print images of variables, show values from diference between variables
* and print image with diference between variables
*
* Denis Eiras
* 21/6/2018 - creation - script compareVariables3.gs
* 13/4/2023 - adapt to WSM5

function main(args)

   if (args='')
      say 'requires arguments: dir1/file.ctl dir2/file.ctl filebasename'
      'quit'
   else
      file1=subwrd(args,1)
      file2=subwrd(args,2)
      filebasename=subwrd(args,3)
   endif
 c say "running with args " args
   fileout=compare(file1,file2,filebasename)

'quit'

function compare(file1,file2,filebasename)
  'reinit'
  openFile(file1,'ctl')
  'set gxout print'
  'set prnopts %s 5 1'
  'q ctlinfo'
  allVars =''
  difmaxima=1.E-6
  incdenom=1.E-12
  resultaux = result
  say resultaux

  i=6
  varsNum=0
  while(varsNum=0)
    tmp = sublin ( resultaux, i )
    i=i+1
    if(subwrd(tmp,1)='vars')
      varsNum = subwrd(tmp,2)
    endif
  endwhile

  openFile(file2,'ctl')
  'set gxout print'
  'set prnopts %s 5 1'
  lineVars=i-1
    ztRange()
    say 'Dimensoes xmax: '_xmax', ymax: '_ymax', zmax: '_zmax
  while(i<=lineVars+varsNum)
    tmp = sublin(resultaux, i)
    var = subwrd(tmp,1)
    zlevels = subwrd(tmp,2)
    if(zlevels=0)
      zlevels=1
    endif
    equalsStringPos=find(var,"=")
    if(equalsStringPos>0)
      var=substr(var,1,equalsStringPos-1)
    endif
    'set x 1 '_xmax
    'set y 1 '_ymax
    say ''
    z=1
    say '=========================================================================================='
    say 'Variavel a ser comparada ===> 'var', Níveis: = 'zlevels
    say
    'set gxout shaded'
    while(z<=zlevels)
      t=1
      while(t<=_tmax)
        'set z 'z
        'set t 't

* Este trecho comentado serve para mostrar a variação das variáveis com relação ao seu range
*        'c'
*        'd 'var
*        cada=sublin(result, 1)
*        'cbarn'
*        'draw title Variavel comparada 'var
*        'printim 'var'_comparada.png'
*        if(subwrd(cada,1)='Constant')
*          varmin=0
*          varmax=0
*        else
*          varmin=subwrd(cada,2)
*          varmax=subwrd(cada,4)
*        endif
*        say '===> Variavel a ser comparada: Min = 'varmin', Max = 'varmax
*        'c'
*        'd 'varmax'-'varmin
*        cada=sublin(result, 1)
*        variacao=subwrd(cada,4)
*        say 'variacao da variavel a ser comparada = 'variacao
*        say

*        'c'
*        'd 'var'.2'
*        cada=sublin(result, 1)
*        if(subwrd(cada,1)='Data')
*          say 'Variavel nova nao encontrada'
*          break
*        endif
*        'cbarn'
*        'draw title Variavel nova 'var
*        'printim 'var'_nova.png'
*        if(subwrd(cada,1)='Constant')
*          varmin2=0
*          varmax2=0
*        else
*          varmin2=subwrd(cada,2)
*          varmax2=subwrd(cada,4)
*        endif
*        say '===> Variavel nova: Min = 'varmin2', Max = 'varmax2
*        'c'
*        'd 'varmax2'-'varmin2
*        cada=sublin(result, 1)
*        variacao2=subwrd(cada,4)
*        say 'variacao da variavel nova = 'variacao2
*        say

        'set gxout shaded'
        'c'

*        d (rain.1-rrp.2)/(rain.1+1.E-12)

        'd abs(('var'-'var'.2)/('var'+'incdenom'))'
        cada=sublin(result, 1)
        'cbarn'
        'draw title Diferenca das variaveis 'var
        if(subwrd(cada,1)='Constant')
          say 'Nao ha diferencas'
        else
          varmax=subwrd(cada,4)
          if(varmax>difmaxima)
            say '===> ATENÇÃO !!! Diferença na variavel 'var' em z='z', t='t' foi de 'varmax' e superou 'difmaxima
           'printim 'var'_z'z'_t't'__diff___'filebasename'.png'

           'c'
           'd 'var
           'set gxout contour'
           'd 'var'.2'
           'cbarn'
           'draw title Contorno 'var'.2 = arq. 2' 
           'printim 'var'_z'z'_t't'__contourn___'filebasename'.png'
          endif
        endif
        t=t+1
      endwhile
      z=z+1
    endwhile
    i= i + 1
  endwhile
  'close 2'
  'close 1'
return fileout

function openFile(fileToOpen,fileExt)
  if(fileExt='.nc')
    'sdfopen 'fileToOpen
  else
    'open 'fileToOpen
  endif
  msg('Arquivo 'fileToOpen' aberto')
return

function msg(messg)
  say '>>>>> 'messg
return

function ztRange()
  'q file'
  tmp = sublin ( result, 5 )
  _xmax = subwrd(tmp,3)
  _ymax = subwrd(tmp,6)
  _zmin = 1
  _zmax = subwrd(tmp,9)
  _tmin = 1
  _tmax = subwrd(tmp,12)
return

function find( str, char )
  ntmp = math_strlen( str )
  i = 1
  while( i <= ntmp )
    tmp = substr(str,i,1)
    if( tmp = char )
      return i
    endif
    i = i + 1
  endwhile
return -1




