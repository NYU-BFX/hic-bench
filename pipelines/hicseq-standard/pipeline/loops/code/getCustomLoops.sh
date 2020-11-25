#!/bin/bash

inpdir=$1
min_anchordist=$2
max_anchordist=$3
min_activity=$4
min_qvalue=$5
inputLoops=$6
topLoops=$7
outdir=$inpdir

################################################### EXAMPLE ### ####################################################################################################################################
## ./code/getCustomLoops.sh results/loops.by_sample.res_10kb_md20kb_qval01/filter.by_sample.mapq_20/align.by_sample.bowtie2/CD34_2190-k27ac 20000 10000000 0 0.1 loops_filtered_nobias_cpm 110000 ##
####################################################################################################################################################################################################

#inpdir=/gpfs/home/rodrij92/LOUCY-HiChIP-aug25/hic-bench/pipelines/hicseq-standard/pipeline/loops/results/loops.by_sample.res_10kb_md20kb_qval01/filter.by_sample.mapq_20/align.by_sample.bowtie2/CD34_2190-k27ac
#min_anchordist=20000
#max_anchordist=10000000
#min_activity=0
#min_qvalue=0.1
#inputLoops=loops_filtered_nobias_cpm
#topLoops=110000
#outdir=$inpdir

awk 'NR>1' ${outdir}/${inputLoops}.tsv | cut -f 7 > ${outdir}/qval.txt

paste ${outdir}/${inputLoops}.bedpe ${outdir}/qval.txt > ${outdir}/loops_labeled_qval.bedpe

awk -v m=${min_anchordist} -v M=${max_anchordist} -v c=${min_activity} -v mqval=${min_qvalue} '($5-$2)>=m && ($5-$2)<=M && $7>=c && $8 <= mqval' ${outdir}/loops_labeled_qval.bedpe | cut -f 1-7 > ${outdir}/loops_labeled.bedpe

mv ${outdir}/loops_filtered_nobias_cpm.bedpe  ${outdir}/loops_filtered_nobias_cpm_original.bedpe 
sort -k7,7 -r ${outdir}/loops_labeled.bedpe | head -n $topLoops > ${outdir}/loops_filtered_nobias_cpm.bedpe
rm -f ${outdir}/qval.txt

echo "Done."
