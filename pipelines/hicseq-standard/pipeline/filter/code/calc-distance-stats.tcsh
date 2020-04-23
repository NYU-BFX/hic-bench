#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: calc-distance-stats FILTERED-REG-GZ-FILE
##

if ( ("$1" == "--help") || ($#argv != 1) ) then
  grep '^##' $0
  exit
endif

set reg = $1

set maxd = 2500000

scripts-send2err "Calculating read pair distance statistics..."

echo "Chromosome,Distance(kb),Count"
cat $reg | gunzip | awk '$2==$6' | gtools-regions strand -v -op + | gtools-regions sort | cut -f2 | awk '{print $1"\t"$0}' | gtools-regions dist | awk -v d=$maxd '$2<=d' | awk '{printf "%s %ld\n",$1,$2/1000}' | sort | uniq -c | awk '{print $2","$3","$1}' 

