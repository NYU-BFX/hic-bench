#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: hicseq-loops-fithic.tcsh OUTPUT-DIR PARAM-SCRIPT HIC-REG-FILES GENOME BRANCH OBJECTS
##

if ($#argv != 6) then
  grep '^##' $0
  exit
endif

set outdir = $1
set params = $2
set reg = ($3)      # *.reg.gz files
set genome = $4     # hg19
set branch = $5
set objects = ($6)

# if objects is empty, use all objects in the branch
if ("$objects" == "") set objects = `cd $branch; ls -1d *`

# get ds-accepted-intra read-pairs count (will be used for contactCount CPM normalization)
set matrix_filtered_branch = `ls -d ../matrix-filtered/results/matrix-filtered.by_sample.res_*/*/* | awk '{if(NR<2) print $1}'` # I have to revise this because: if more than 1 aligner params is used in the hicbench-run, here we will only use stats.tsv of the 1st one in the list.

@ intra_reads = 0
foreach obj ($objects)
  @ x = `cat $matrix_filtered_branch/$obj/stats.tsv | grep ^ds-accepted-intra | cut -f2`
  @ intra_reads = $intra_reads + $x
end
echo "ds-accepted-intra-reads = $intra_reads"

# Run the parameter script
source $params

# Set up genome information
set genome_dir = inputs/genomes/$genome
set genome_file = $genome_dir/genome/bowtie2.index/genome.fa.fai
set chromosomes = `cat $genome_dir/genome.bed | cut -f1 | grep -vw "$chrom_excluded"`
set n_chromosomes = $#chromosomes 
@ n_chromosomes --                                           # this is because we want to run an array job --array=0:(n-1)

# Create uncompressed version of mapped read pairs
echo "Uncompressing filtered reads..." | scripts-send2err
cat $reg | gunzip >! $outdir/filtered.reg

# Enter object's directory
set main_dir = `echo ${cwd}`

# Create bin files
echo "Creating bin files..." | scripts-send2err
mkdir -p $outdir/bins
cd $outdir/bins
windowMaker -g $main_dir/$genome_file -w $winsize >! x.bed
awk '{printf "%s\t%d\t%d\t%s_%.0f\n",$1,$2+1,$3,$1,($2+$3)*0.5}' x.bed >! k.bed
awk '{print>$1}' k.bed
rm -f x.bed
cd $main_dir

# Convert filtered.reg bed file to bedpe format
echo "Converting reg to bedpe format..." | scripts-send2err
mkdir -p $outdir/bedpe
awk -v OFS='\t' '{if ($2 == $6) {print $2,$4,$5,$6,$8,$9,$1,".",$3,$7}}' $outdir/filtered.reg >! $outdir/bedpe/intra_converted.reg
cd $outdir/bedpe/
awk '{print >$1}' intra_converted.reg #splits the bedpe file by chromosome
cd $main_dir

# Call loops for each chromosome in a separate folder
set job_dir = $outdir/__jdata
mkdir $job_dir
set jid = `sbatch --array=0-$n_chromosomes --output="$job_dir/job.%a.out" --error="$job_dir/job.%a.err" ./code/scripts-loops-fithic-chr.sh $outdir $outdir/bins $winsize "$chromosomes"`
set jid = `echo $jid | sed 's/.* //'`
echo $jid >! $job_dir/job.id
scripts-send2err "Waiting for job array [$jid] to complete..."
scripts-qsub-wait "$jid"

# Concatenate chromosome loops in one file
scripts-send2err "Combining chromosome loops into one file..."
cd $outdir 

# unfiltered loops
cat chr*/loops_unfiltered.tsv >! temp.tsv
awk 'NR <= 1 || \!/fragment/' temp.tsv >! all_loops_unfiltered.tsv
gzip -f all_loops_unfiltered.tsv
rm -f temp.tsv

# filtered loops
cat chr*/loops_filtered.tsv >! temp.tsv
awk 'NR <= 1 || \!/fragment/' temp.tsv >! all_loops_filtered.tsv
awk -v var="$winsize" '{ if ((NR>1)) print $1"\t"($2-var/2)"\t"($2+var/2)"\t"$3"\t"($4-var/2)"\t"($4+var/2)"\t"$5}' all_loops_filtered.tsv >! all_loops_filtered.bedpe
rm -f temp.tsv

# create IGV junction format (loops-like)
awk '{if(NR>1) print $1"\t"$2"\t"$4"\t\.\t1\.0"}' all_loops_filtered.tsv | sed -e '1itrack graphType=junctions' | sort -k2 -n >! all_loops_filtered.igv.bed

# Create CPM normalized loops files
awk -v var="$intra_reads" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5/(var/1000000)"\t"$6"\t"$7"\t"$8"\t"$9}' all_loops_filtered.tsv >! all_loops_filtered_cpm.tsv
awk -v var="$intra_reads" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7/(var/1000000)}' all_loops_filtered.bedpe >! all_loops_filtered_cpm.bedpe

# Clean up
rm -fr bedpe bins chr* slurm* filtered.reg 
