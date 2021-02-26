#!/bin/bash
#SBATCH -J chipEmul
#SBATCH --mem=20G
#SBATCH --time=3:00:00
#SBATCH -N 1

inpdir=$1
species=$2

module load macs2/2.1.1 
module load bedtools/2.27.1
module load ucscutils/374

#### HOW TO RUN  ############################################################################
# chipEmul.sh <filtered.reg directory> <genome>   # genome options: hs (human) or mm (mice) #
#############################################################################################

#### RUN EXAMPLE ####################################################################################################
# sbatch ./code/chipEmul.sh results/filter.by_sample.mapq_20_mindist0/align.by_sample.bowtie2/CUTLL1-DMSO-k27ac/ hs #
#####################################################################################################################

genome=`grep -w "^genome" $inpdir/job.vars.tsv | cut -f2`
genome_dir=`grep -w "^genome_dir" $inpdir/job.vars.tsv | cut -f2`
idx=${genome_dir}/bowtie2.index/${genome}.fa.fai

echo $inpdir
echo $genome
echo $genome_dir
echo $idx

# gen short paired reads with insert size < 1kb
gunzip -c $inpdir/filtered.reg.gz > $inpdir/filtered.reg
sed  's/ /\t/g' $inpdir/filtered.reg > $inpdir/x
awk '{print $2"\t"$4"\t"$4+49"\t"$6"\t"$8"\t"$8+49"\t"".""\t"".""\t"$3"\t"$7}' $inpdir/x > $inpdir/x1  ###I add 49 if it is a 50bp read
awk '{ if ($1 == $4) { print } }' $inpdir/x1 > $inpdir/x2

awk -v OFS='\t' '{		# fix positions so start < end
  	if ($2 < $5)
		print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10;
	else
		print $1,$5,$6,$4,$2,$3,$7,$8,$9,$10;
}' $inpdir/x2 > $inpdir/x2_fix
mv $inpdir/x2_fix $inpdir/x2

awk '{ if($5-$3 <= 1000) { print }}' $inpdir/x2 > $inpdir/x3
nReads=`cat $inpdir/x3 | wc -l`
echo "total short paired end reads: "$nReads

sclFactor=`echo "scale=10 ; 1000000 / $nReads" | bc`
echo "scaling factor: "$sclFactor

# prepare file for genomeCov
cat $inpdir/x3 > $inpdir/shortpaired.bed
cut -f1,2,6 $inpdir/shortpaired.bed > $inpdir/1shortpaired1.bed 
sort -k1,1 -k2,2n $inpdir/1shortpaired1.bed > $inpdir/sort_1shortpaired1.bed
fgrep -v "chrM" $inpdir/sort_1shortpaired1.bed > $inpdir/temp.bed		# to avoid an indexing error
mv $inpdir/temp.bed $inpdir/sort_1shortpaired1.bed

# make bedgraph
genomeCoverageBed -bg -scale $sclFactor -i $inpdir/sort_1shortpaired1.bed -g $idx > $inpdir/hichip_chip_scaled.bedGraph

# make bigwig
bedGraphToBigWig $inpdir/hichip_chip_scaled.bedGraph $idx $inpdir/hichip_chip_scaled.bw

# call peaks
macs2 callpeak -f BEDPE --keep-dup all -g $species -t $inpdir/x3 --outdir $inpdir -n hichip_chip --broad
macs2 callpeak -f BEDPE --keep-dup all -g $species -t $inpdir/x3 --outdir $inpdir -n hichip_chip

# clean up
rm -f $inpdir/x $inpdir/x1 $inpdir/x2 $inpdir/x3 $inpdir/filtered.reg $inpdir/hichip_chip.bedGraph $inpdir/*shortpaired*
