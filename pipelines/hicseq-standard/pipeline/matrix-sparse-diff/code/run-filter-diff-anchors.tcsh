#!/bin/tcsh

##
## USAGE: run-filter-diff-anchors.tcsh RESULTS-DIR
##

if ( ("$1" == "--help") || ($#argv != 1) ) then
  grep '^##' $0
  exit
endif

set results = $1

module load r

set X = `find $results -name diff-anchors.csv | sed 's/.[^/]\+$//' | sort`

foreach out ($X)
  echo $out | sed 's/.*\///'
  Rscript ./code/filter-diff-anchors.r --min-dist=10000 --min-val=4.0 $out/diff-anchors.csv > ! $out/filtered-diff-anchors.csv
  #wc -l $out/*diff-anchors.csv
end


