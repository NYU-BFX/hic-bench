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
set tss_path = `readlink -f $TSS_genes`

set main_dir = `echo ${cwd}`
cd $outdir
set sampleName = `basename "$PWD"`

## Reorder columns (filtered.reg -> filtered.bed = HiCsummary format from Homer) ##
awk '{print $1 "\t" $2 "\t" $4 "\t" $3 "\t" $6 "\t" $9 "\t" $7}' filtered.reg >! filtered.temp
awk '{ if ($2 == $5) { print } }' filtered.temp >! filtered.bed

##### RUN HOMER PIPELINE #####
## Create TAG directory ##
makeTagDirectory TAG -format HiCsummary filtered.bed

## Principal Component Analysis of Hi-C Data ##
runHiCpca.pl pca_HKgenesFix TAG -res ${resolution} -pc 2 -cpu 8 -active $hkgene_path
runHiCpca.pl pca_tssFix TAG -res ${resolution} -cpu 8 -active $tss_path

## Fix sign ##
intersectBed -a pca_HKgenesFix.PC1.bedGraph -b $hkgene_path -c > pca_HKgenesFix.PC1_counts.bed
intersectBed -a pca_tssFix.PC1.bedGraph -b $tss_path -c > pca_tssFix.PC1_counts.bed
Rscript $main_dir/code/scripts-compartments-metrics.r pca_HKgenesFix.PC1_counts.bed pca_tssFix.PC1_counts.bed "${sampleName}"

## Find HiC compartments ##
findHiCCompartments.pl pca_HKgenesFix.PC1.txt >! pca_HKgenesFix_Acompartments.txt
findHiCCompartments.pl pca_HKgenesFix.PC1.txt -opp >! pca_HKgenesFix_Bcompartments.txt

awk '{print $2"\t"$3"\t"$4}' pca_HKgenesFix_Acompartments.txt >! pca_HKgenesFix_Acompartments.bed
awk '{print $2"\t"$3"\t"$4}' pca_HKgenesFix_Bcompartments.txt >! pca_HKgenesFix_Bcompartments.bed

# Clean up
rm -fr filtered.temp filtered.reg filtered.bed TAG pca_HKgenesFix_Acompartments.txt pca_HKgenesFix_Bcompartments.txt pca_tssFix.PC1* pca_HKgenesFix.PC1_counts.bed
