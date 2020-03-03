#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: hicseq-tracks-washu.tcsh OUTPUT-DIR HIC-REG-FILES GENOME-BED BIN-SIZE
##
## Example: ./hicseq-tracks-washu.tcsh mESC_J1-HindIII "mESC_J1-HindIII-rep1*/filtered_reads.reg+" genomes/mm10/genome.bed 50000 
##

if ($#argv != 4) then
  grep '^##' $0
  exit
endif

set out = $1
set r = ($2)
set g = $3     # genomes/mm10/genome.bed
set b = $4     # 50000

# Create genome matrix
cat $r | gunzip | gtools-hic bin -v --bin-size $b -g $g --matrix | gtools-hic convert >! $out/track.washu.tsv

# Compress and index
bgzip $out/track.washu.tsv
tabix -p bed $out/track.washu.tsv.gz



