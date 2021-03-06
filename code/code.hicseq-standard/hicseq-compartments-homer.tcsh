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

#shuf -n 5000000 $outdir/filtered.reg >! $outdir/filtered2.reg  # only for testing
#rm -f $outdir/filtered.reg					# only for testing
#mv $outdir/filtered2.reg $outdir/filtered.reg			# only for testing

# Enter object's directory
set hkgene_path = `readlink -f $HK_genes`
set tss_path = `readlink -f $TSS_genes`

if ($active_mark != FALSE) then
	echo $active_mark 
	set hkgene_path = `readlink -f $active_mark`
endif

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
runHiCpca.pl pca_tssFix TAG -res ${resolution} -pc 2 -cpu 8 -active $tss_path

mv pca_HKgenesFix.PC1.bedGraph pca_HKgenesFix.PC1.PC2.bedGraph
mv pca_tssFix.PC1.bedGraph pca_tssFix.PC1.PC2.bedGraph

sed -n '/PC2/q;p' pca_HKgenesFix.PC1.PC2.bedGraph > pca_HKgenesFix.PC1.bedGraph
sed -n -e '/PC2/,$p' pca_HKgenesFix.PC1.PC2.bedGraph > pca_HKgenesFix.PC2.bedGraph
sed -n '/PC2/q;p' pca_tssFix.PC1.PC2.bedGraph > pca_tssFix.PC1.bedGraph

## Fix sign ##
intersectBed -a pca_HKgenesFix.PC1.bedGraph -b $hkgene_path -c > pca_HKgenesFix.PC1_counts.bed
intersectBed -a pca_tssFix.PC1.bedGraph -b $tss_path -c > pca_tssFix.PC1_counts.bed
Rscript $main_dir/code/scripts-compartments-metrics.r pca_HKgenesFix.PC1_counts.bed pca_tssFix.PC1_counts.bed "${sampleName}"

## Find HiC compartments ##
mv pca_HKgenesFix.PC1.txt pca_HKgenesFix.PC1.PC2.txt
cut -f1-6 pca_HKgenesFix.PC1.PC2.txt  > pca_HKgenesFix.PC1.txt
findHiCCompartments.pl pca_HKgenesFix.PC1.txt >! pca_HKgenesFix_Acompartments.txt
findHiCCompartments.pl pca_HKgenesFix.PC1.txt -opp >! pca_HKgenesFix_Bcompartments.txt

awk '{print $2"\t"$3"\t"$4}' pca_HKgenesFix_Acompartments.txt >! A_compartments.bed
awk '{print $2"\t"$3"\t"$4}' pca_HKgenesFix_Bcompartments.txt >! B_compartments.bed

mv pca_HKgenesFix.PC1.bedGraph compartments.scores.bedGraph 

## Clean up
rm -fr filtered.temp filtered.reg filtered.bed TAG pca_HKgenesFix_Acompartments.txt pca_HKgenesFix_Bcompartments.txt pca_tssFix.PC1* pca_HKgenesFix.PC1_counts.bed pca_HKgenesFix.PC1.PC2.bedGraph pca_HKgenesFix.PC1.txt
#rm -f pca_HK* metrics*

fgrep -v "Error" job.err > temp.txt; mv temp.txt job.err # remove homer's error messages that doesn't seem to be relevant 
