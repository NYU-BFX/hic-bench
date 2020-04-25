#!/bin/tcsh

##
## USAGE: run-fix
##

if ( ("$1" == "--help") || ($#argv != 0) ) then
  grep '^##' $0
  exit
endif

module load r

set X = results-BALL-new/matrix-sparse-diff.by_sample.win_10kb/matrix-sparse.by_sample.unit_100bp.maxd_5Mb/filter.by_sample.mapq_20/align.by_sample.bowtie2/BALL-*/*

foreach out ($X)
  echo $out | sed 's/.*\///'
  Rscript ./code/filter-diff-anchors.r --min-dist=10000 --min-val=4.0 $out/diff-anchors.csv > ! $out/filtered-diff-anchors.csv
  wc -l $out/*diff-anchors.csv
end


