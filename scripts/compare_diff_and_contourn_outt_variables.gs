'reinit'

'open ../dataout/ref_g.ctl'
'open ../dataout__compare_local_modifications/ref_g.ctl'

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
'printim outt'var'_diff.png'

say ""
say "Comparing shaded (original) x contourn (actual)"
say ""
'c'
'set gxout shaded'
'd outt'var'.1'
'set gxout contour'
'd outt'var'.2'
'printim outt'var'_contourn.png'

var=var+1
endwhile

'quit'

