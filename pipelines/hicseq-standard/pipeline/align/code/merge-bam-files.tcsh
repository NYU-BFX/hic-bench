#!/bin/tcsh

##
## USAGE: merge-bam-files.tcsh R1-BAM R2-BAM GENOME-FA
##

if ( ("$1" == "--help") || ($#argv != 3) ) then
  grep '^##' $0
  exit
endif

set R1 = $1
set R2 = $2
set fa = $3

samtools view $R1 >! $R1.sam
samtools view $R2 >! $R2.sam

paste -d'\n' $R1.sam $R2.sam | samtools view -T $fa -b1 - 

rm -f $R1.sam $R2.sam

