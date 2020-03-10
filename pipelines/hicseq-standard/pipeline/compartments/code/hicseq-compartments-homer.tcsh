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

# Convert to bed format
awk ' BEGIN { OFS="\t"; strand["-"]="1"; strand["+"]="0" } {                    \
    if ($2 < $6)                                                                \
        print $1, strand[$3], $2, $4, 0, strand[$7], $6, $8, 1, 0, 1;           \
    else                                                                        \
        print $1, strand[$7], $6, $8, 1, strand[$3], $2, $4, 0, 0, 1;           \
}' $outdir/filtered.reg | sort -k3,3d -k7,7d >! $outdir/filtered.bed

##### RUN HOMER PIPELINE #####
## Create TAG directory ##
makeTagDirectory $outdir/TAG -format HiCsummary $outdir/filtered.bed

## Principal Component Analysis of Hi-C Data ##
runHiCpca.pl $outdir/pca_activeMarkFix $outdir/TAG -res ${res} -cpu 8 -active ${HK_genes}

## Find HiC compartments ##
findHiCCompartments.pl $outdir/pca_activeMarkFix.PC1.txt >! $outdir/pca_activeMarkFix_Acompartments.txt
findHiCCompartments.pl $outdir/pca_activeMarkFix -opp >! $outdir/pca_activeMarkFix_Bcompartments.txt

awk '{print $2"\t"$3"\t"$4}' $outdir/pca_activeMarkFix_Acompartments.txt >! $outdir/pca_activeMarkFix_Acompartments.bed
awk '{print $2"\t"$3"\t"$4}' $outdir/pca_activeMarkFix_Bcompartments.txt >! $outdir/pca_activeMarkFix_Vcompartments.txt

