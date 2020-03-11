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

ls -l $HK_genes

# Create uncompressed version of mapped read pairs
cat $reg | gunzip >! $outdir/filtered.reg

## Reorder columns (filtered.reg -> filtered.bed = HiCsummary format from Homer) ##
awk '{print $1 "\t" $2 "\t" $4 "\t" $3 "\t" $6 "\t" $9 "\t" $7}' $outdir/filtered.reg > $outdir/filtered.temp
awk '{ if ($2 == $5) { print } }'  $outdir/filtered.temp > $outdir/filtered.bed

##### RUN HOMER PIPELINE #####
## Create TAG directory ##
makeTagDirectory $outdir/TAG -format HiCsummary $outdir/filtered.bed

## Principal Component Analysis of Hi-C Data ##
runHiCpca.pl $outdir/pca_activeMarkFix $outdir/TAG -res ${resolution} -cpu 8 -active ${HK_genes}

## Find HiC compartments ##
findHiCCompartments.pl $outdir/pca_activeMarkFix.PC1.txt >! $outdir/pca_activeMarkFix_Acompartments.txt
findHiCCompartments.pl $outdir/pca_activeMarkFix.PC1.txt -opp >! $outdir/pca_activeMarkFix_Bcompartments.txt

awk '{print $2"\t"$3"\t"$4}' $outdir/pca_activeMarkFix_Acompartments.txt >! $outdir/pca_activeMarkFix_Acompartments.bed
awk '{print $2"\t"$3"\t"$4}' $outdir/pca_activeMarkFix_Bcompartments.txt >! $outdir/pca_activeMarkFix_Bcompartments.bed

rm -fr $outdir/filtered.temp $outdir/filtered.reg $outdir/filtered.bed $outdir/TAG $outdir/pca_activeMarkFix_Acompartments.txt $outdir/pca_activeMarkFix_Bcompartments.txt
