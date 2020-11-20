#!/bin/bash

LOOPS=$1
K27AC=$2
OUTDIR=$3
OUTNAME=$4

awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$1":"$2":"$3,$4":"$5":"$6}' ${OUTDIR}/${LOOPS} > ${OUTDIR}/loops_temp.txt

cut -f 1-3,8 ${OUTDIR}/loops_temp.txt > ${OUTDIR}/a1.txt
cut -f 4-6,9 ${OUTDIR}/loops_temp.txt > ${OUTDIR}/a2.txt

bedtools intersect -a ${OUTDIR}/a1.txt -b ${K27AC} -wao | cut -f 1-4,8 | sort -k1,1 -k2,2n | bedtools merge -i stdin -c 5,5,5,4 -o mean,sum,max,distinct | sort -k7,7b > ${OUTDIR}/a1_k27ac_signal.tsv

bedtools intersect -a ${OUTDIR}/a2.txt -b ${K27AC} -wao | cut -f 1-4,8 | sort -k1,1 -k2,2n | bedtools merge -i stdin -c 5,5,5,4 -o mean,sum,max,distinct | sort -k7,7b > ${OUTDIR}/a2_k27ac_signal.tsv

sort -k8,8b ${OUTDIR}/loops_temp.txt | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$1":"$2":"$3":"$4":"$5":"$6}'  > ${OUTDIR}/loops_a1_sorted.txt
sort -k9,9b ${OUTDIR}/loops_temp.txt | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$1":"$2":"$3":"$4":"$5":"$6}' > ${OUTDIR}/loops_a2_sorted.txt

echo "chr1 start1 end1 chr2 start2 end2 cpm a1_mean.enh.act a1_sum.enh.act a1_max.enh.act a2_mean.enh.act a2_sum.enh.act a2_max.enh.act sum.loop.enh.act mean.loop.enh.act" | tr ' ' '\t' >> ${OUTDIR}/${OUTNAME}
join -1 8 -2 7 ${OUTDIR}/loops_a1_sorted.txt ${OUTDIR}/a1_k27ac_signal.tsv | sort -k9,9b > ${OUTDIR}/a1_data.tsv
join -1 9 -2 7 ${OUTDIR}/a1_data.tsv ${OUTDIR}/a2_k27ac_signal.tsv | sort -u -k10,10b | tr ' ' '\t' | cut -f 3-9,14-16,20-22 | sed 's/\t\./\t0/g' | awk -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,($9+$12),(($9+$12)/2)}' >> ${OUTDIR}/${OUTNAME}

# clean up
rm -f ${OUTDIR}/loops_temp.txt ${OUTDIR}/a1.txt ${OUTDIR}/a2.txt ${OUTDIR}/a1_data.tsv ${OUTDIR}/a1_k27ac_signal.tsv ${OUTDIR}/a2_k27ac_signal.tsv
