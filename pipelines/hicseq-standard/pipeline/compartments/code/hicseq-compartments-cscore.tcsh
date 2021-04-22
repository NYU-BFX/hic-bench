#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: hicseq-compartments-cscore.tcsh OUTPUT-DIR PARAM-SCRIPT HIC-REG-FILES GENOME
##

if ($#argv != 6) then
  grep '^##' $0
  exit
endif

set outdir = $1
set params = $2
set reg = ($3)      # *.reg.gz files
set inpdirs = ($4) # inpdirs/filter/results/f*/a*/*, length is equal to the number of samples used
set objects = ($5) # type-1-rep1 type-1-rep2, length is equal to the number of samples used
set genome = $6     # hg19

# Run the parameter script
source $params

# Create uncompressed version of mapped read pairs
cat $reg | gunzip >! $outdir/filtered.reg
# cat $reg | gunzip | head -n 100000 >! $outdir/filtered.reg

# Enter object's directory
set hkgene_path = `readlink -f $HK_genes`
set tss_path = `readlink -f $TSS_genes`

if ($active_mark != FALSE) then
	echo $active_mark 
        set hkgene_path = `readlink -f $active_mark`
endif

# echo $hkgene_path
# echo $tss_path

# set main_dir = `echo ${cwd}`
# cd $outdir
set main_dir = $cwd
set sampleName = `basename "$outdir"`

## Reorder columns (filtered.reg -> filtered.bed = HiCsummary format from Homer) ##
awk '{print $1 "\t" $2 "\t" $4 "\t" $3 "\t" $6 "\t" $9 "\t" $7}' $outdir/filtered.reg >! $outdir/filtered.temp
awk '{ if ($2 == $5) { print } }' $outdir/filtered.temp >! $outdir/filtered.bed

##### RUN CSCORE PIPELINE #####

# This compartment step is *independent* of matrix-filtered. No need to reference bin_size. Instead, use $resolution (default 100000).

module unload bedtools
module load bedtools
bedtools makewindows -g ${main_dir}/inputs/genomes/$genome/bowtie2.index/$genome.fa.fai -w $resolution > $outdir/$genome.tiled.$resolution.temp.bed
module unload bedtools
# User might wonder why we have to manually add +1 to start column.
# The primary reason is sloppy code in cscoretool: in order for the final bedgraph output to conform to 0 based coordinate.
awk '{$2 = $2 + 1 ; print $1"\t"$2"\t"$3}' $outdir/$genome.tiled.$resolution.temp.bed > $outdir/$genome.tiled.$resolution.bed


# cat `/gpfs/data/abl/home/choh09/programs/CscoreTool/CscoreTool1.1 --help`

# echo "Dollar inpdirs is : $inpdirs"
# echo "Length of inpdirs is: $#inpdirs"
# echo "Dollar objects is : $objects"
# echo "Length of objects is: $#objects"
# echo $outdir
# echo $mindist

# run CscoreTool in parallel
# Example run:
# /gpfs/data/abl/home/choh09/programs/CscoreTool/CscoreTool1.1 mm10.tiled.100k.20200604.bed intermediate.2.20200604.bed ./chr1/c1 10 1000000 chr1

# Create temporary directory, for each chromosome:
mkdir $cwd/$outdir/__jdata
set jid = 
set chrs = `awk '{print $1}' ${main_dir}/inputs/genomes/$genome/bowtie2.index/$genome.fa.fai`
# echo $chrs
# echo $#chrs

foreach chr ($chrs)
  scripts-send2err "Processing $chr..."
  
  # estimate required memory
  set mem = "32G"
  scripts-send2err "requested memory = $mem"

  # run ScoreTool inside
  set jpref = $outdir/each.`echo $chr | sed 's/\.[^.]\+$//'`
  # echo $jpref
  # scripts-create-path $jpref
  mkdir $cwd/$jpref
  # echo "whgusdn $chr" > $jpref/$chr.txt
  set jid = ($jid `scripts-qsub-run $jpref 1 $mem $cwd/code/CscoreTool1.1 $cwd/$outdir/$genome.tiled.$resolution.bed $cwd/$outdir/filtered.bed $cwd/$outdir/each.${chr}/${chr} 10 $mindist $chr`)

end

# wait until all jobs are completed
scripts-send2err "Waiting until all jobs are completed..."
scripts-qsub-wait "$jid"

# Sign flipping must occur *after* all bedgraph files are made and the processes are completed.
scripts-send2err "Because of availability of housekeeping genes, do not expect the sign change on chrs M, X, Y. Take extreme care, if you need results on them."
scripts-send2err "Sign flipping step here."
foreach chr ($chrs)
  ## Fix sign ##
  # Invoke R code here!!
  Rscript $main_dir/code/scripts-compartments-cscore-flip-signs.R $outdir/each.${chr}/${chr}_cscore.bedgraph $hkgene_path
  ## Fix sign ## end
end

# Aggregate all *cscore.bedgraph files.
cat `ls $outdir/*/*_cscore.bedgraph` > $outdir/CscoreTool.scores.temp.bedGraph

# Need only one "track type bedGraph" at the top. Only once.
echo 'track type="bedGraph"' > $outdir/compartments.scores.bedGraph
grep -v "^track" $outdir/CscoreTool.scores.temp.bedGraph >> $outdir/compartments.scores.bedGraph

# track type="bedGraph" name="/gpfs/data/abl/home/choh09/Work/05_May/Effie_visualize/pipelines/hicseq-standard/pipeline/compartments/results/compartments.by_sample.cscore.res_100kb/filter.by_sample.mapq_20/align.by_sample.bowtie2/MEF-untreated-Arima-rep2/__jdata/each.chr1/chr1_cscore.bedgraph"

## Find HiC compartments ##
# findHiCCompartments.pl pca_HKgenesFix.PC1.txt >! pca_HKgenesFix_Acompartments.txt
# findHiCCompartments.pl pca_HKgenesFix.PC1.txt -opp >! pca_HKgenesFix_Bcompartments.txt

# awk '{print $2"\t"$3"\t"$4}' pca_HKgenesFix_Acompartments.txt >! pca_HKgenesFix_Acompartments.bed
# awk '{print $2"\t"$3"\t"$4}' pca_HKgenesFix_Bcompartments.txt >! pca_HKgenesFix_Bcompartments.bed

# # Clean up
rm -rf $outdir/filtered.temp $outdir/filtered.reg $outdir/filtered.bed $outdir/$genome.tiled.$resolution.temp.bed $outdir/$genome.tiled.$resolution.bed $outdir/$genome.tiled.$resolution.bed $outdir/CscoreTool.scores.temp.bedGraph
rm -fr $outdir/each.chr*
