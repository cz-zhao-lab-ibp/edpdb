#!/bin/tcsh

# this script is for large scale screening of 3D motif. 
# written by X.C. Zhang, on 090524

source ~/edpdb/edpdb.csh 

cd /hosts/xtal/zhangc/pdb/pdb/

cat <<$eof > find_clique.edp
init; atom zn ; if (-i status <= 0 ) goto finish
group zn
init; ca ; more ; swap ; group nprt
! exclude proteins of multiple-metals active sites 
init; atom zn fe mg co cd mn cu na k | load nprt ; doit zn s 1 0 5.0 ; if ( -i status >= 1 ) goto finish 

init; residue Cys His Asp Glu | side ; exclude atom c* ; group sc
init; residue His | atom nd1 ce1 cd2 ne2 ; group His

parameter cm = z
step1:  
! do not let chain z/j mess up
init; chain \$(cm) ; if (-i status <= 0 ) goto step2
if ( -s cm == j ) goto finish
parameter cm = j  
goto step1

step2: 
read template_pdb.tmp \$(cm)
@find_zn 
init; load his1 ; exclude load xx4 ; doit ; if (-i status <= 1) goto finish
{ residue Asp Glu | atom o* | load sc }
group haz
init; chain \$(cm) | atom oe2 ne2 zn ; swap haz 

setenv verbose 2 ; clique haz 3 1.2 1.2 20 1 1 ; if ( -i status <= 0 ) goto finish
file 
finish:
quit
$eof

cat <<$eof1 > find_zn.edp
init; atom xx1 xx2 xx3 xzn ; group xx1
init; atom xx1 xx3 xx0 xzn ; group xx2 
init; atom xx4 ; group xx4
init; atom xx* xzn ; group xall
init; group his1 
init; atom zn ; group zn0 

loop:
parameter zn1 zn0 id exit
init; zone \$(zn1) ; group zn1 
doit sc s 1 0 2.6 
init; load scr ; more | ca ; doit ; if (-i status ^= 3 ) goto loop
init; load zn1 
doit sc s 1 1.0 2.2 ; if (-i status == 3 ) goto ok2
doit sc s 1 1.0 2.3 ; if (-i status == 3 ) goto ok2
doit sc s 1 1.0 2.4 ; if (-i status == 3 ) goto ok2
doit sc s 1 1.0 2.5 ; if (-i status == 3 ) goto ok2
doit sc s 1 1.0 2.6 ; if (-i status ^= 3 ) goto loop
ok2:
init; load scr ; group lgd3
init; load xx1 ; over lgd3 ; load xall ; rtn file 
init; atom xzn ; doit zn1 s 1 0 0.5 ; if (-i status == 1 ) goto ok1
init; load xx2 ; over lgd3 ; load xall ; rtn file
init; atom xzn ; doit zn1 s 1 0 0.5 ; if (-i status ^= 1 ) goto loop
ok1:
init; load his ; doit xx4 s 1 4.0 8.0 zn1 s 1 2.0 3.0 a 90 140 ; if (-i status <= 0 ) goto loop
init; load scr his1 ; group his1
goto loop
$eof1

cat <<$eof2 > template_pdb.tmp 
ATOM   1229  OE2 GLU A 156     138.409 254.705 124.058  1.00 76.18
ATOM   1797  ND1 HIS A 225     138.029 254.466 126.574  1.00 75.83
ATOM   1799  NE2 HIS A 225     137.393 254.772 128.599  1.00 74.77
ATOM   3243 ZN    ZN A 501     132.038 253.386 132.573  1.00 65.05
ATOM   2251  XX1 xxx A   2      35.654  83.154 132.474  1.00 78.65
ATOM   2624  XX2 xxx A   2      33.834  85.909 130.755  1.00 71.98
ATOM   2631  XX3 xxx A   2      33.500  85.941 134.137  1.00 78.28
ATOM   2624  XX0 xxx A   2      33.834  85.909 130.755  1.00 71.98
ATOM   6474  XX4 xxx A   2      36.932  86.593 132.964  1.00 69.17
ATOM   3245  XZN xxx A   2      34.852  85.307 132.759  1.00 75.56
END
$eof2

#edpdb 1oht.pdb -i find_clique 
#exit 

foreach i ( ????.pdb  ) 
	edpdb $i -i find_clique
end

rm find_*.edp  *_.txt template_pdb.tmp
exit

