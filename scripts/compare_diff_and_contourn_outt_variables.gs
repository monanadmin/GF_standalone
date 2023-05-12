
function main(args)

   if (args='')
      say 'requires arguments: testingdir/file.ctl referencedir/file.ctl'
      'quit'
   else
      file1=subwrd(args,1)
      file2=subwrd(args,2)
   endif
 c say "running with args " args
   fileout=compare(file1,file2)

'quit'

function compare(file1,file2)
'reinit'

'open 'file1
'open 'file2

var=1
while(var<=3)

say "Comparing out"var" variable"
say "-------------------------"
say ""
say "Diff"
say ""
'c'
'set gxout shaded'
'd outt'var'.2 - outt'var'.1'
*'cbar'
'printim 'file1'___outt'var'_diff'.png'

say ""
say "Comparing shaded (original) x contourn (actual)"
say ""
'c'
'set gxout shaded'
'd outt'var'.1'
'set gxout contour'
'd outt'var'.2'
'printim 'file1'___outt'var'_contourn.png'

var=var+1
endwhile

'quit'

