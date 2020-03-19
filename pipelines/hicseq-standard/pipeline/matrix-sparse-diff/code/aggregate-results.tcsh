#!/bin/tcsh

##
## USAGE: aggregate-results.tcsh
##

if ($#argv != 0) then
  grep '^##' $0
  exit
endif

set ppp = `pwd`
set out = `readlink -f .`/BALL-stats.csv
cd BALL
set PATIENTS = `ls -1d *`
head -1 $PATIENTS[1]/*/stats.csv | sed 's/^/Patient,/' >! $out
foreach p ($PATIENTS)
  echo "Processing patient $p..."
  set s = $p/*/stats.csv
  cat $s | skipn 1 | awk -F',' '$8>=100 && ($7-$6>100) && ($7>300)' | sed "s/^/$p,+/" >> $out
  cat $s | skipn 1 | awk -F',' '$8>=100 && ($6-$7>100) && ($6>300)' | sed "s/^/$p,-/" >> $out
end

cd $ppp
cat BALL-stats.csv | cut -d',' -f1,2 | sort -u | cut -d',' -f2 | sort | uniq -cd | awk '$1>=3' | sort -n

