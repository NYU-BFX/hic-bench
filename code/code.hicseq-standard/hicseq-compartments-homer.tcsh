#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: hicseq-compartments-homer.tcsh OUTPUT-DIR PARAM-SCRIPT HIC-REG-FILES GENOME
##

if ($#argv != 4) then
  grep '^##' $0
  exit
endif

set outdir = $1
set params = $2
set reg = ($3)      # *.reg.gz files
set genome = $4     # hg19

# Run the parameter script
source $params

# Create uncompressed version of mapped read pairs
cat $reg | gunzip >! $outdir/filtered.reg

# Enter object's directory
set hkgene_path = `readlink -f $HK_genes`
cd $outdir

## Reorder columns (filtered.reg -> filtered.bed = HiCsummary format from Homer) ##
awk '{print $1 "\t" $2 "\t" $4 "\t" $3 "\t" $6 "\t" $9 "\t" $7}' filtered.reg >! filtered.temp
awk '{ if ($2 == $5) { print } }' filtered.temp >! filtered.bed

##### RUN HOMER PIPELINE #####
## Create TAG directory ##
makeTagDirectory TAG -format HiCsummary filtered.bed

## Principal Component Analysis of Hi-C Data ##
runHiCpca.pl pca_activeMarkFix TAG -res ${resolution} -cpu 8 -active $hkgene_path

## Find HiC compartments ##
findHiCCompartments.pl pca_activeMarkFix.PC1.txt >! pca_activeMarkFix_Acompartments.txt
findHiCCompartments.pl pca_activeMarkFix.PC1.txt -opp >! pca_activeMarkFix_Bcompartments.txt

awk '{print $2"\t"$3"\t"$4}' pca_activeMarkFix_Acompartments.txt >! pca_activeMarkFix_Acompartments.bed
awk '{print $2"\t"$3"\t"$4}' pca_activeMarkFix_Bcompartments.txt >! pca_activeMarkFix_Bcompartments.bed

# Clean up
rm -fr filtered.temp filtered.reg filtered.bed TAG pca_activeMarkFix_Acompartments.txt pca_activeMarkFix_Bcompartments.txt


