#!/bin/tcsh 
# you may specify this shell program as the default viewer of PDB file 
#  

cd `dirname $1`
echo "PROGRAM `alias edpdb`"
echo "FILE   "  $1 
echo "DIR     `pwd`"
echo 

edpdb `basename $1` -i edpini.edp

if( $status >< 0) then 

cd $HOME
echo "PROGRAM `alias edpdb`"
echo "FILE   "  $1 
echo "DIR     `pwd`"
echo 

edpdb $1 -i edpini.edp
endif

echo ; echo "type <enter> to close this window" ; line
exit 0
