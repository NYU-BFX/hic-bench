#!/bin/bash
#SBATCH -J peakCall_2
#SBATCH --mem=80G
#SBATCH --time=24:00:00
#SBATCH -N 1

#### HOW TO RUN  #####################################################
# sbatch --array=1-n scripts-peakCall.sh  <branch> ; n = num samples #
######################################################################

branch=$1

sample_idx=${SLURM_ARRAY_TASK_ID}
echo $sample_idx
echo $branch
sampleName=$(ls -l $branch | fgrep -v "total" |  awk "NR==${sample_idx}" | awk '{print $9}')
outdir=${branch}/${sampleName}
echo $outdir

# Detect species
genome_prefix=`cat $outdir/job.vars.tsv | fgrep "genome" | awk 'NR==1' | cut -f2 | cut -c1-2`

if [ $genome_prefix = 'hg' ] 
then
	species='hs'
	echo 'Detected human (hs) species.'
else
	species='mm'
	echo 'Detected mouse (mm) species.'
fi	

module load bedtools/2.27.1
module load ucscutils/374

# Get fasta index file
genome=`grep -w "^genome" $outdir/job.vars.tsv | cut -f2`
genome_dir=`grep -w "^genome_dir" $outdir/job.vars.tsv | cut -f2`
idx=${genome_dir}/bowtie2.index/${genome}.fa.fai

# gen short paired reads with insert size < 1kb
gunzip -c $outdir/filtered.reg.gz > $outdir/filtered.reg
sed  's/ /\t/g' $outdir/filtered.reg > $outdir/x
awk '{print $1,"\t"$2"\t"$4"\t"$4+49"\t"$6"\t"$8"\t"$8+49"\t"".""\t"".""\t"$3"\t"$7}' $outdir/x > $outdir/x1  ###I add 49 if it is a 50bp read
awk '{ if ($2 == $5) { print } }' $outdir/x1 > $outdir/x2

awk -v OFS='\t' '{		# fix positions so start < end
  	if ($3 < $6)
		print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11;
	else
		print $1,$5,$6,$4,$2,$3,$7,$8,$9,$10,$11;
}' $outdir/x2 > $outdir/x2_fix
mv $outdir/x2_fix $outdir/x2

awk '{ if($6-$4 <= 1000) { print }}' $outdir/x2 > $outdir/short_readPairs.tsv
cut -f1 $outdir/short_readPairs.tsv > $outdir/chip_read_ids.txt

nReads=`cat $outdir/short_readPairs.tsv | wc -l`
echo "total short paired end reads: "$nReads

# make chip.bam
module unload python
module unload macs2/2.1.1
module load python/gcc/3.6.5
module load samtools

#align_branch=`echo $outdir | cut -d"/" -f3-4`
align_branch=`echo $outdir | awk -F '/filter/' '{print $2}'| cut -d '/' -f3-4`
align_dir=inpdirs/align/results/$align_branch

samtools merge $outdir/allReads.bam $align_dir/bowtie_results/bwt2/*/*.bwt2pairs.bam
samtools sort $outdir/allReads.bam -o $outdir/sorted.bam
samtools index $outdir/sorted.bam
python ./code/extract_reads.py -b $outdir/sorted.bam -n $outdir/chip_read_ids.txt -o $outdir/extracted.bam
samtools view -H -o $outdir/chip.sam $outdir/extracted.bam
samtools view -o $outdir/extracted.sam $outdir/extracted.bam
awk '!seen[$1]++' $outdir/extracted.sam >> $outdir/chip.sam
samtools view -S -b $outdir/chip.sam > $outdir/chip.bam
samtools sort $outdir/chip.bam -o $outdir/chip_sorted.bam
samtools index $outdir/chip_sorted.bam

# compute scaling factor for cpm normalization
sclFactor=`echo "scale=10 ; 1000000 / $nReads" | bc`
echo "scaling factor: "$sclFactor

# prepare file for genomeCov
cut -f 2-11 $outdir/short_readPairs.tsv > $outdir/shortpaired.bed
cut -f1,2,6 $outdir/shortpaired.bed > $outdir/1shortpaired1.bed 
sort -k1,1 -k2,2n $outdir/1shortpaired1.bed > $outdir/sort_1shortpaired1.bed
fgrep -v "chrM" $outdir/sort_1shortpaired1.bed > $outdir/temp.bed		# to avoid an indexing error
mv $outdir/temp.bed $outdir/sort_1shortpaired1.bed

# make bedgraph
genomeCoverageBed -bg -scale $sclFactor -i $outdir/sort_1shortpaired1.bed -g $idx > $outdir/chip_scaled.bedGraph

# make bigwig
bedGraphToBigWig $outdir/chip_scaled.bedGraph $idx $outdir/chip_scaled.bw

# call peaks
module unload python
module load macs2/2.1.1

macs2 callpeak -f BEDPE --keep-dup all -g $species -t $outdir/shortpaired.bed --outdir $outdir -n chip --broad
#macs2 callpeak -f BEDPE --keep-dup all -g $species -t $outdir/shortpaired.bed --outdir $outdir -n chip

# clean up
rm -f $outdir/x $outdir/x1 $outdir/x2 $outdir/filtered.reg $outdir/chip.bedGraph $outdir/*shortpaired* $outdir/extracted.* $outdir/chip.sam $outdir/chip_peaks.gappedPeak $outdir/chip_peaks.xls $outdir/sorted.bam* $outdir/chip_summits.bed
rm -f $outdir/short_readPairs.tsv $outdir/chip_scaled.bedGraph $outdir/chip.bam 
