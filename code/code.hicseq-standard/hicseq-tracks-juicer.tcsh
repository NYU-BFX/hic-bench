#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: hicseq-tracks-juicer.tcsh OUTPUT-DIR HIC-REG-FILES GENOME
##
## Example: ./hicseq-tracks-washu.tcsh mESC_J1-HindIII "mESC_J1-HindIII-rep1*/filtered_reads.reg+" mm10 
##

if ($#argv != 3) then
  grep '^##' $0
  exit
endif

set outdir = $1
set reg = ($2)      # *.reg.gz files
set genome = $3     # mm10

# Create uncompressed version of mapped read pairs
cat $reg | gunzip >! $outdir/filtered.reg

# Convert to bed format
awk ' BEGIN { OFS="\t"; strand["-"]="1"; strand["+"]="0" } {
    if ($2 < $6)
        print $1, strand[$3], $2, $4, 0, strand[$7], $6, $8, 1, 0, 1;
    else
        print $1, strand[$7], $6, $8, 1, strand[$3], $2, $4, 0, 0, 1;
}' $outdir/filtered.reg | sort -k3,3d -k7,7d >! $outdir/filtered.bed

# Generate .hic file
java -Xms512m -Xmx20480m -jar juicer_tools.jar pre $outdir/filtered.bed $outdir/filtered.hic "$genome"

# Cleanup
rm -f $outdir/filtered.reg $outdir/filtered.bed

