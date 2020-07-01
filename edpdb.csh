#To run EDPDB, one needs to include the following into one's 	
# .cshrc file.							
								
# variables for EDPDB						
setenv  EDPDBIN    /xtal/edpdb				
setenv  EDP_DATA   $EDPDBIN/data					

alias	edpdb      $EDPDBIN/edpdb_v13a

exit

#for .bashrc, use the following

declare -x EDPDBIN=/usr/local/edpdb
declare -x EDP_DATA=$EDPDBIN/data
alias edpdb=$EDPDBIN/edpdb_v13a
