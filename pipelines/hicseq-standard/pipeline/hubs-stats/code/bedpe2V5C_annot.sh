#!/bin/bash
BEDPE=$1
K27AC=$2
TSS=$3
OUTDIR=$4

module load bedtools/2.27.1

rm -fr ${OUTDIR}
mkdir -p ${OUTDIR}

awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$1":"$2":"$3":"$4":"$5":"$6}' $BEDPE > ${OUTDIR}/loops_labeled.bedpe
cut -f 1-3,7-8 ${OUTDIR}/loops_labeled.bedpe | awk -v OFS="\t" '{print $1,$2,$3,$1":"$2":"$3,$4,$5}' > ${OUTDIR}/a1.bed 
cut -f 4-6,7-8 ${OUTDIR}/loops_labeled.bedpe | awk -v OFS="\t" '{print $1,$2,$3,$1":"$2":"$3,$4,$5}' > ${OUTDIR}/a2.bed 

bedtools intersect -a ${K27AC} -b ${OUTDIR}/a1.bed -wo | sort -k10b,10 > ${OUTDIR}/k27ac_a1_intersect.txt
bedtools intersect -a ${K27AC} -b ${OUTDIR}/a2.bed -wo | sort -k10b,10 > ${OUTDIR}/k27ac_a2_intersect.txt
bedtools intersect -a ${TSS} -b ${OUTDIR}/a1.bed -wo | sort -k12b,12 > ${OUTDIR}/tss_a1_intersect.txt
bedtools intersect -a ${TSS} -b ${OUTDIR}/a2.bed -wo | sort -k12b,12 > ${OUTDIR}/tss_a2_intersect.txt

join -j10 ${OUTDIR}/k27ac_a1_intersect.txt ${OUTDIR}/k27ac_a2_intersect.txt | awk '$9!=$19'> ${OUTDIR}/EE_loops.txt
join -j12 ${OUTDIR}/tss_a1_intersect.txt ${OUTDIR}/tss_a2_intersect.txt | awk '$11!=$23' > ${OUTDIR}/PP_loops.txt
join -1 10 -2 12 ${OUTDIR}/k27ac_a1_intersect.txt ${OUTDIR}/tss_a2_intersect.txt | awk '$9!=$21' > ${OUTDIR}/EP_loops.txt
join -1 12 -2 10 ${OUTDIR}/tss_a1_intersect.txt ${OUTDIR}/k27ac_a2_intersect.txt | awk '$11!=$21' > ${OUTDIR}/PE_loops.txt

echo "Total number of EE loops:"
wc -l ${OUTDIR}/EE_loops.txt
echo "Total number of PP loops:"
wc -l ${OUTDIR}/PP_loops.txt
echo "Total number of EP loops:"
wc -l ${OUTDIR}/EP_loops.txt
echo "Total number of PE loops:"
wc -l ${OUTDIR}/PE_loops.txt

# get bedpe
awk -v OFS="\t" '{print $6,$7,$8,$16,$17,$18,$10,$5,$15}' ${OUTDIR}/EE_loops.txt > ${OUTDIR}/EE_loops.bedpe
awk -v OFS="\t" '{print $8,$9,$10,$20,$21,$22,$12,$5,$17}' ${OUTDIR}/PP_loops.txt > ${OUTDIR}/PP_loops.bedpe
awk -v OFS="\t" '{print $6,$7,$8,$18,$19,$20,$10,$5,$15}' ${OUTDIR}/EP_loops.txt > ${OUTDIR}/EP_loops.bedpe
awk -v OFS="\t" '{print $8,$9,$10,$18,$19,$20,$12,$5,$17}' ${OUTDIR}/PE_loops.txt > ${OUTDIR}/PE_loops.bedpe
cat ${OUTDIR}/EP_loops.bedpe ${OUTDIR}/PE_loops.bedpe > ${OUTDIR}/EP_PE_loops.bedpe
cat ${OUTDIR}/EP_PE_loops.bedpe ${OUTDIR}/PP_loops.bedpe > ${OUTDIR}/EP_PE_PP_loops.bedpe
cat ${OUTDIR}/EE_loops.bedpe ${OUTDIR}/EP_PE_PP_loops.bedpe > ${OUTDIR}/all_loops.bedpe

awk -v OFS="\t" '{print $1":"$2":"$3":"$4":"$5":"$6,$1,$2,$3,$4,$5,$6,$7,$8,$9}' ${OUTDIR}/EE_loops.bedpe | sort -u -k1,1 | cut -f 2-8 > ${OUTDIR}/EE_loops_uniq.bedpe 
awk -v OFS="\t"	'{print $1":"$2":"$3":"$4":"$5":"$6,$1,$2,$3,$4,$5,$6,$7,$8,$9}' ${OUTDIR}/PP_loops.bedpe | sort -u -k1,1 | cut -f 2-8 > ${OUTDIR}/PP_loops_uniq.bedpe
awk -v OFS="\t" '{print $1":"$2":"$3":"$4":"$5":"$6,$1,$2,$3,$4,$5,$6,$7,$8,$9}' ${OUTDIR}/EP_loops.bedpe | sort -u -k1,1 | cut -f 2-8 > ${OUTDIR}/EP_loops_uniq.bedpe
awk -v OFS="\t" '{print $1":"$2":"$3":"$4":"$5":"$6,$1,$2,$3,$4,$5,$6,$7,$8,$9}' ${OUTDIR}/PE_loops.bedpe | sort -u -k1,1 | cut -f 2-8 > ${OUTDIR}/PE_loops_uniq.bedpe
awk -v OFS="\t" '{print $1":"$2":"$3":"$4":"$5":"$6,$1,$2,$3,$4,$5,$6,$7,$8,$9}' ${OUTDIR}/EP_PE_loops.bedpe | sort -u -k1,1 | cut -f 2-8 > ${OUTDIR}/EP_PE_loops_uniq.bedpe
awk -v OFS="\t" '{print $1":"$2":"$3":"$4":"$5":"$6,$1,$2,$3,$4,$5,$6,$7,$8,$9}' ${OUTDIR}/EP_PE_PP_loops.bedpe | sort -u -k1,1 | cut -f 2-8 > ${OUTDIR}/EP_PE_PP_loops_uniq.bedpe
awk -v OFS="\t" '{print $1":"$2":"$3":"$4":"$5":"$6,$1,$2,$3,$4,$5,$6,$7,$8,$9}' ${OUTDIR}/all_loops.bedpe | sort -u -k1,1 | cut -f 2-8 > ${OUTDIR}/all_loops_uniq.bedpe

echo "Total number of unique loops with EE, PP, EP or PE annotation:"
wc -l ${OUTDIR}/all_loops_uniq.bedpe

echo "Total number of unique loops with PP, EP or PE annotation:"
wc -l ${OUTDIR}/EP_PE_PP_loops_uniq.bedpe

# get loops in virtual5C.csv format
awk -v OFS="\t" '{print $8":"$9,$8,$9,$1,($2+$3)/2,($5+$6)/2,$5-$2,"-",$7,"-","-","0","0","-","-"}' ${OUTDIR}/all_loops.bedpe | sort -u -k1,1 > ${OUTDIR}/all_loops_v5cFormat.tsv

# add reversed duplicated loops (for v5cFormat consistency)
awk -v OFS="\t" '{print $3":"$2,$3,$2,$4,$6,$5,$7,$8,$9,$10,$11,$12,$13,$14,$15}' ${OUTDIR}/all_loops_v5cFormat.tsv | sort -u -k1,1 > ${OUTDIR}/all_loops_rev_v5cFormat.tsv

# get final v5cFormat loops file
echo "Source.anchor.label,Target.anchor.label,Chromosome,Source.anchor.position,Target.anchor.position,Anchor.distance,Count,CPK2B,VP.total.count,Count.sum,pvalue,fdr,p.expected,p.observed" >> ${OUTDIR}/all_loops_wRev_v5cFormat.csv
cat ${OUTDIR}/all_loops_v5cFormat.tsv ${OUTDIR}/all_loops_rev_v5cFormat.tsv | sort -u -k1,1 | cut -f 2-15 | tr '\t' ',' >> ${OUTDIR}/all_loops_wRev_v5cFormat.csv

# clean up
rm -f ${OUTDIR}/a1.bed ${OUTDIR}/a2.bed ${OUTDIR}/*_intersect.txt ${OUTDIR}/*_loops.txt ${OUTDIR}/loops_labeled.bedpe

echo "Done."
