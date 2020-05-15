#!/bin/tcsh

##
## USAGE: test-symmetry.tcsh INPDIR
##

if ( ("$1" == "--help") || ($#argv != 1) ) then
  grep '^##' $0
  exit
endif

set inpdir = $1
set out = `new_temp`  # out.tsv

set t = `new_temp`

set inpfiles = $inpdir/diff-anchors.csv
if (! -e $inpfiles) set inpfiles = $inpdir/chr*/diff-anchors.csv
cat $inpfiles | awk -F, '$6!=0' | cut -d',' -f1,2,6-8 | sed 's/,-/,/' | tr ',' '\t' >! $t
( cat $t ; cat $t | cols -t 1 0 2 3 4 ) | sort | grep -v label | awk '$1<$2' >! $out 

head $out
set n = `cat $out | sort -u | wc -l`
set k = `cat $out | sort | uniq -c | awk '$1!=2' | wc -l`
echo "$k/$n"
echo "$k/$n" | bc -l

rm -f $t $out

