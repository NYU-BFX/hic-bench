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

# Count intra-reads
set intra_reads = `cat $outdir/filtered.reg | awk '$2 == $6' | wc -l`
echo "ds-accepted-intra-reads = $intra_reads" | scripts-send2err

# Enter object's directory
set main_dir = `echo ${cwd}`

# Create bin files
echo "Creating bin files..." | scripts-send2err
mkdir -p $outdir/bins
cd $outdir/bins
module unload bedtools
module load bedtools
bedtools makewindows -g $main_dir/$genome_file -w $winsize >! x.bed
module unload bedtools
awk '{printf "%s\t%d\t%d\t%s_%.0f\n",$1,$2+1,$3,$1,($2+$3)*0.5}' x.bed >! k.bed
awk '{print>$1}' k.bed
rm -f x.bed
cd $main_dir

# Convert filtered.reg bed file to bedpe format
echo "Converting reg to bedpe format..." | scripts-send2err
mkdir -p $outdir/bedpe

awk -v OFS='\t' '{                                                        \
  	if ($2 == $6 && $5 < $8)     			                  \
		print $2,$4,$5,$6,$8,$9,$1,".",$3,$7;                     \
	else if ($2 == $6 && $5 > $8)                                     \
		print $2,$8,$9,$6,$4,$5,$1,".",$3,$7;                     \
}' $outdir/filtered.reg >! $outdir/bedpe/intra_converted.reg

cd $outdir/bedpe/
awk '{print >$1}' intra_converted.reg #splits the bedpe file by chromosome
cd $main_dir

### Identify loops with FitHic ###
# Call loops for each chromosome in a separate folder
echo "Calling loops for each chromosome..." | scripts-send2err
set job_dir = $outdir/__jdata
mkdir -p $job_dir
set jid = `sbatch --array=0-$n_chromosomes --output="$job_dir/job.%a.out" --error="$job_dir/job.%a.err" ./code/scripts-loops-fithic-chr.sh $outdir $outdir/bins $winsize "$chromosomes" "$qval"`
set jid = `echo $jid | sed 's/.* //'`
echo $jid >! $job_dir/job.id
scripts-send2err "Waiting for job array [$jid] to complete..."
scripts-qsub-wait "$jid"

# Concatenate chromosome loops in one file
scripts-send2err "Combining chromosome loops into one file..."
cd $outdir 

# unfiltered loops raw (bias)
cat chr*/loops_unfiltered_bias_raw.tsv >! temp.tsv
awk 'NR <= 1 || \!/fragment/' temp.tsv >! loops_unfiltered_bias_raw.tsv
awk -v var="$winsize" '{ if ((NR>1)) print $1"\t"($2+1-var/2)"\t"($2+var/2)"\t"$3"\t"($4+1-var/2)"\t"($4+var/2)"\t"$5}' loops_unfiltered_bias_raw.tsv >! loops_unfiltered_bias_raw.bedpe
rm -f temp.tsv

# unfiltered loops (no bias)
cat chr*/loops_unfiltered_nobias_raw.tsv >! temp.tsv
awk 'NR <= 1 || \!/fragment/' temp.tsv >! loops_unfiltered_nobias_raw.tsv
awk -v var="$winsize" '{ if ((NR>1)) print $1"\t"($2+1-var/2)"\t"($2+var/2)"\t"$3"\t"($4+1-var/2)"\t"($4+var/2)"\t"$5}' loops_unfiltered_nobias_raw.tsv >! loops_unfiltered_nobias_raw.bedpe
rm -f temp.tsv

# Create CPM normalized unfiltered loops files (bias)
awk -v var="$intra_reads" '{                                                              \
	if ( NR == 1 )                                                                    \
		print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9                  \
	else                                                                              \
		print $1"\t"$2"\t"$3"\t"$4"\t"$5/(var/1000000)"\t"$6"\t"$7"\t"$8"\t"$9    \
}' loops_unfiltered_bias_raw.tsv >! loops_unfiltered_bias_cpm.tsv
awk -v var="$winsize" '{ if ((NR>1)) print $1"\t"($2+1-var/2)"\t"($2+var/2)"\t"$3"\t"($4+1-var/2)"\t"($4+var/2)"\t"$5}' loops_unfiltered_bias_cpm.tsv >! loops_unfiltered_bias_cpm.bedpe

# Create CPM normalized unfiltered loops files (no bias)
awk -v var="$intra_reads" '{                                                              \
        if ( NR == 1 )                                                                    \
                print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9                  \
        else                                                                              \
                print $1"\t"$2"\t"$3"\t"$4"\t"$5/(var/1000000)"\t"$6"\t"$7"\t"$8"\t"$9    \
}' loops_unfiltered_nobias_raw.tsv >! loops_unfiltered_nobias_cpm.tsv
awk -v var="$winsize" '{ if ((NR>1)) print $1"\t"($2+1-var/2)"\t"($2+var/2)"\t"$3"\t"($4+1-var/2)"\t"($4+var/2)"\t"$5}' loops_unfiltered_nobias_cpm.tsv >! loops_unfiltered_nobias_cpm.bedpe


# filtered loops (bias)
cat chr*/loops_filtered_bias_raw.tsv >! temp.tsv
awk -v min="$mindist" -v max="$maxdist" 'NR == 1 || \!/fragment/ && ($4+1-$2) < max && ($4+1-$2) > min' temp.tsv >! loops_filtered_bias_raw.tsv
awk -v var="$winsize" '{ if ((NR>1)) print $1"\t"($2+1-var/2)"\t"($2+var/2)"\t"$3"\t"($4+1-var/2)"\t"($4+var/2)"\t"$5}' loops_filtered_bias_raw.tsv >! loops_filtered_bias_raw.bedpe
rm -f temp.tsv

# filtered loops (no bias)
cat chr*/loops_filtered_nobias_raw.tsv >! temp.tsv
awk -v min="$mindist" -v max="$maxdist" 'NR == 1 || \!/fragment/ && ($4+1-$2) < max && ($4+1-$2) > min' temp.tsv >! loops_filtered_nobias_raw.tsv
awk -v var="$winsize" '{ if ((NR>1)) print $1"\t"($2+1-var/2)"\t"($2+var/2)"\t"$3"\t"($4+1-var/2)"\t"($4+var/2)"\t"$5}' loops_filtered_nobias_raw.tsv >! loops_filtered_nobias_raw.bedpe
rm -f temp.tsv

# create IGV junction format (loops-like)
#awk '{if(NR>1) print $1"\t"$2"\t"$4"\t\.\t1\.0"}' loops_filtered_bias_raw.tsv | sed -e '1itrack graphType=junctions' | sort -k2 -n >! loops_filtered_bias.igv.bed
#awk '{if(NR>1) print $1"\t"$2"\t"$4"\t\.\t1\.0"}' loops_filtered_nobias_raw.tsv | sed -e '1itrack graphType=junctions' | sort -k2 -n >! loops_filtered_nobias.igv.bed

# Create CPM normalized filtered loops files (bias)
awk -v var="$intra_reads" '{                                                              \
	if ( NR == 1 )                                                                    \
		print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9                  \
	else                                                                              \
		print $1"\t"$2"\t"$3"\t"$4"\t"$5/(var/1000000)"\t"$6"\t"$7"\t"$8"\t"$9    \
}' loops_filtered_bias_raw.tsv >! loops_filtered_bias_cpm.tsv

awk -v var="$intra_reads" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7/(var/1000000)}' loops_filtered_bias_raw.bedpe >! loops_filtered_bias_cpm.bedpe

# Create CPM normalized filtered loops files (no bias)
awk -v var="$intra_reads" '{                                                              \
        if ( NR == 1 )                                                                    \
                print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9                  \
        else                                                                              \
                print $1"\t"$2"\t"$3"\t"$4"\t"$5/(var/1000000)"\t"$6"\t"$7"\t"$8"\t"$9    \
}' loops_filtered_nobias_raw.tsv >! loops_filtered_nobias_cpm.tsv

awk -v var="$intra_reads" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7/(var/1000000)}' loops_filtered_nobias_raw.bedpe >! loops_filtered_nobias_cpm.bedpe


### scale the cpm values so the mean cpm value is 1 ###

# no bias
set mean_cpm = `awk 'NR>1' loops_filtered_nobias_cpm.bedpe | cut -f7 | awk '{ total += $1 } END { print total/NR }'`
set mfc = `echo $mean_cpm | awk '{print 1/$1}'`
echo 'The original mean cpm value of the nobias loops is = '$mean_cpm
echo 'The correction factor is = '$mfc

awk -v OFS="\t" -v f="$mfc" '{printf "%s\t%d\t%d\t%s\t%d\t%d\t%.4f\n", $1,$2,$3,$4,$5,$6,$7*f}' loops_filtered_nobias_cpm.bedpe >! loops_filtered_nobias_cpmScaled.bedpe
set mean_scaled = `awk 'NR>1' loops_filtered_nobias_cpmScaled.bedpe | cut -f7 | awk '{ total += $1 } END { print total/NR }'`
echo 'The mean cpm value of the nobias loops after correction is = '$mean_scaled

# bias
set mean_cpm = `awk 'NR>1' loops_filtered_bias_cpm.bedpe | cut -f7 | awk '{ total += $1 } END { print total/NR }'`
set mfc = `echo $mean_cpm | awk '{print 1/$1}'`
echo 'The original mean cpm value of the bias loops is = '$mean_cpm
echo 'The correction factor is = '$mfc

awk -v OFS="\t" -v f="$mfc" '{printf "%s\t%d\t%d\t%s\t%d\t%d\t%.4f\n", $1,$2,$3,$4,$5,$6,$7*f}' loops_filtered_bias_cpm.bedpe >! loops_filtered_bias_cpmScaled.bedpe
set mean_scaled = `awk 'NR>1' loops_filtered_bias_cpmScaled.bedpe | cut -f7 | awk '{ total += $1 } END { print total/NR }'`
echo 'The mean cpm value of the bias loops after correction is = '$mean_scaled


###  select top loops and perform cpm scaling on the top loops selected ###
echo $top_loops
if ($top_loops != "FALSE") then
	# no bias
	sort -nr -k7,7 loops_filtered_nobias_cpm.bedpe | head -n $top_loops >! loops_filtered_nobias_cpm_top${top_loops}.bedpe
	set mean_cpm = `awk 'NR>1' loops_filtered_nobias_cpm_top${top_loops}.bedpe | cut -f7 | awk '{ total += $1 } END { print total/NR }'`
	set mfc = `echo $mean_cpm | awk '{print 1/$1}'`
	echo 'The original mean cpm value of the nobias top loops is = '$mean_cpm
	echo 'The correction factor is = '$mfc

	awk -v OFS="\t" -v f="$mfc" '{printf "%s\t%d\t%d\t%s\t%d\t%d\t%.4f\n", $1,$2,$3,$4,$5,$6,$7*f}' loops_filtered_nobias_cpm_top${top_loops}.bedpe >! loops_filtered_nobias_top${top_loops}_cpmScaled.bedpe
	set mean_scaled = `awk 'NR>1' loops_filtered_nobias_top${top_loops}_cpmScaled.bedpe | cut -f7 | awk '{ total += $1 } END { print total/NR }'`
	echo 'The mean cpm value of the nobias top loops after correction is = '$mean_scaled

	# bias
	sort -nr -k7,7 loops_filtered_bias_cpm.bedpe | head -n $top_loops >! loops_filtered_bias_cpm_top${top_loops}.bedpe
	set mean_cpm = `awk 'NR>1' loops_filtered_bias_cpm_top${top_loops}.bedpe | cut -f7 | awk '{ total += $1 } END { print total/NR }'`
	set mfc = `echo $mean_cpm | awk '{print 1/$1}'`
	echo 'The original mean cpm value of the bias top loops is = '$mean_cpm
	echo 'The correction factor is = '$mfc

	awk -v OFS="\t" -v f="$mfc" '{printf "%s\t%d\t%d\t%s\t%d\t%d\t%.4f\n", $1,$2,$3,$4,$5,$6,$7*f}' loops_filtered_bias_cpm_top${top_loops}.bedpe >! loops_filtered_bias_top${top_loops}_cpmScaled.bedpe
	set mean_scaled = `awk 'NR>1' loops_filtered_bias_top${top_loops}_cpmScaled.bedpe | cut -f7 | awk '{ total += $1 } END { print total/NR }'`
	echo 'The mean cpm value of the bias top loops after correction is = '$mean_scaled
endif

### Create QC plots ###

# Get random unfiltered 500k loops
head -n 1 loops_unfiltered_bias_raw.tsv > loops_unfiltered_bias_raw_shuf500k.tsv
tail -n +2 loops_unfiltered_bias_raw.tsv | shuf -n 500000 >> loops_unfiltered_bias_raw_shuf500k.tsv

# compress big files
gzip -f loops_unfiltered_bias_cpm.tsv
gzip -f loops_unfiltered_nobias_cpm.tsv
gzip -f loops_unfiltered_bias_raw.tsv
gzip -f loops_unfiltered_nobias_raw.tsv
gzip -f loops_unfiltered_nobias_raw.bedpe
gzip -f loops_unfiltered_nobias_cpm.bedpe
gzip -f loops_unfiltered_bias_raw.bedpe
gzip -f loops_unfiltered_bias_cpm.bedpe

#rm -f loops_unfiltered_bias_raw.tsv loops_unfiltered_nobias_raw.tsv loops_unfiltered_bias_cpm.tsv loops_unfiltered_nobias_cpm.tsv

mkdir -p QC_plots
cd $main_dir
@ n = $n_chromosomes - 1
Rscript ./code/scripts-loops-QC.r $outdir $n

# Clean up
cd $outdir
rm -fr loops_unfiltered_bias_raw_shuf500k.tsv bedpe bins chr* slurm* filtered.reg
