#!/bin/bash
BEDPE=$1
K27AC=$2
TSS=$3
ATAC=$4
accessible_only=$5
tss_extension=$6
promoter_k27ac_only=$7
OUTDIR=$8

#BEDPE=/Users/javrodher/Work/RStudio-PRJs/ESC-MEF-Alex/data/RUN-TEST/esc_fithic-loops-10kb_nobias_cpm_md10.bedpe
#K27AC=/Users/javrodher/Work/RStudio-PRJs/ESC-MEF-Alex/data/RUN-TEST/k27ac_labeled_esc.bed
#TSS=/Users/javrodher/Work/RStudio-PRJs/ESC-MEF-Alex/data/RUN-TEST/3tss.bed
#OUTDIR=/Users/javrodher/Work/RStudio-PRJs/ESC-MEF-Alex/data/hubs-stats-sep24-md10/ES-untreated/bedpe2V5C/test

module load bedtools/2.27.1

rm -fr ${OUTDIR}
mkdir -p ${OUTDIR}

awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$1":"$2":"$3":"$4":"$5":"$6}' $BEDPE > ${OUTDIR}/loops_labeled.bedpe
cut -f 1-3,7-8 ${OUTDIR}/loops_labeled.bedpe | awk -v OFS="\t" '{print $1,$2,$3,$1":"$2":"$3,$4,$5}' > ${OUTDIR}/a1.bed 
cut -f 4-6,7-8 ${OUTDIR}/loops_labeled.bedpe | awk -v OFS="\t" '{print $1,$2,$3,$1":"$2":"$3,$4,$5}' > ${OUTDIR}/a2.bed 

# clean k27ac file
cut -f 1-3,5 ${K27AC} > ${OUTDIR}/k27ac_clean.bed
K27AC=${OUTDIR}/k27ac_clean.bed 

# filter by accessibility

if [[ $accessible_only = "TRUE" ]]
then
	awk -v OFS="\t" -v ext=$tss_extension '{
	if ($6 == "+")
    		print $1,$2-ext,$3,$4
	else
    		print $1,$2,$3+ext,$4
	}' ${TSS} > ${OUTDIR}/tss_extended.bed

	bedtools intersect -a ${K27AC} -b ${ATAC} -wa -u > ${OUTDIR}/k27ac_flt.bed
	bedtools intersect -a ${OUTDIR}/tss_extended.bed -b ${ATAC} -wa -u | sort -k4,4b > ${OUTDIR}/tss_flt_extended.bed
	sort -k4,4b ${TSS} > ${OUTDIR}/tss_sorted.bed 
	join -1 4 -2 4 ${OUTDIR}/tss_sorted.bed ${OUTDIR}/tss_flt_extended.bed | awk -v OFS="\t" '{print $2,$3,$4,$1,$5,$6}' > ${OUTDIR}/tss_flt.bed
	rm -f ${OUTDIR}/tss_sorted.bed ${OUTDIR}/tss_flt_extended.bed ${OUTDIR}/tss_extended.bed
else
	cp ${K27AC} ${OUTDIR}/k27ac_flt.bed
	cp ${TSS} ${OUTDIR}/tss_flt.bed
fi

# filter promoters without k27ac peaks

if [[ $promoter_k27ac_only = "TRUE" ]]
then
	awk -v OFS="\t" -v ext=$tss_extension '{
	if ($6 == "+")
    		print $1,$2-ext,$3,$4
	else
    		print $1,$2,$3+ext,$4
	}' ${OUTDIR}/tss_flt.bed > ${OUTDIR}/tss_extended.bed

	bedtools intersect -a ${OUTDIR}/tss_extended.bed -b ${K27AC} -wa -u | sort -k4,4b > ${OUTDIR}/tss_flt_extended.bed
	sort -k4,4b ${OUTDIR}/tss_flt.bed > ${OUTDIR}/tss_sorted.bed 
	rm -f ${OUTDIR}/tss_flt.bed
	join -1 4 -2 4 ${OUTDIR}/tss_sorted.bed ${OUTDIR}/tss_flt_extended.bed | awk -v OFS="\t" '{print $2,$3,$4,$1,$5,$6}' > ${OUTDIR}/tss_flt.bed
fi

# intersect k27ac and TSS with anchors
bedtools intersect -a ${OUTDIR}/k27ac_flt.bed -b ${OUTDIR}/a1.bed -wo | sort -k10b,10 > ${OUTDIR}/k27ac_a1_intersect.txt
bedtools intersect -a ${OUTDIR}/k27ac_flt.bed -b ${OUTDIR}/a2.bed -wo | sort -k10b,10 > ${OUTDIR}/k27ac_a2_intersect.txt
bedtools intersect -a ${OUTDIR}/tss_flt.bed -b ${OUTDIR}/a1.bed -wo | sort -k12b,12 > ${OUTDIR}/tss_a1_intersect.txt
bedtools intersect -a ${OUTDIR}/tss_flt.bed -b ${OUTDIR}/a2.bed -wo | sort -k12b,12 > ${OUTDIR}/tss_a2_intersect.txt

wc -l ${OUTDIR}/k27ac_a1_intersect.txt
wc -l ${OUTDIR}/k27ac_a2_intersect.txt

# remove k27ac peaks that are in a TSS anchor
cat ${OUTDIR}/k27ac_a1_intersect.txt | sort -k8b,8 > ${OUTDIR}/K1
cat ${OUTDIR}/tss_a1_intersect.txt | sort -k10b,10 > ${OUTDIR}/T1

cat ${OUTDIR}/k27ac_a2_intersect.txt | sort -k8b,8 > ${OUTDIR}/K2
cat ${OUTDIR}/tss_a2_intersect.txt | sort -k10b,10 > ${OUTDIR}/T2

join -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11  -v 1 -1 8 -2 10 ${OUTDIR}/K1 ${OUTDIR}/T1 | tr ' ' '\t' | sort -k10b,10 > ${OUTDIR}/temp; mv ${OUTDIR}/temp ${OUTDIR}/k27ac_a1_intersect.txt
join -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11  -v 1 -1 8 -2 10 ${OUTDIR}/K2 ${OUTDIR}/T2 | tr ' ' '\t' | sort -k10b,10 > ${OUTDIR}/temp; mv ${OUTDIR}/temp ${OUTDIR}/k27ac_a2_intersect.txt

wc -l ${OUTDIR}/k27ac_a1_intersect.txt
wc -l ${OUTDIR}/k27ac_a2_intersect.txt

# get EE,EP,PE,PP loops
join -j10 ${OUTDIR}/k27ac_a1_intersect.txt ${OUTDIR}/k27ac_a2_intersect.txt | awk '$9!=$19'> ${OUTDIR}/EE_loops.txt
join -j12 ${OUTDIR}/tss_a1_intersect.txt ${OUTDIR}/tss_a2_intersect.txt | awk '$11!=$23' > ${OUTDIR}/PP_loops.txt
join -1 10 -2 12 ${OUTDIR}/k27ac_a1_intersect.txt ${OUTDIR}/tss_a2_intersect.txt | awk '$9!=$21' > ${OUTDIR}/EP_loops.txt
join -1 12 -2 10 ${OUTDIR}/tss_a1_intersect.txt ${OUTDIR}/k27ac_a2_intersect.txt | awk '$11!=$21' > ${OUTDIR}/PE_loops.txt

# get bedpe
awk -v OFS="\t" '{print $6,$7,$8,$16,$17,$18,$10,$5,$15}' ${OUTDIR}/EE_loops.txt > ${OUTDIR}/EE_loops.bedpe
awk -v OFS="\t" '{print $8,$9,$10,$20,$21,$22,$12,$5,$17}' ${OUTDIR}/PP_loops.txt > ${OUTDIR}/PP_loops.bedpe
awk -v OFS="\t" '{print $6,$7,$8,$18,$19,$20,$10,$5,$15}' ${OUTDIR}/EP_loops.txt > ${OUTDIR}/EP_loops.bedpe
awk -v OFS="\t" '{print $8,$9,$10,$18,$19,$20,$12,$5,$17}' ${OUTDIR}/PE_loops.txt > ${OUTDIR}/PE_loops.bedpe
cat ${OUTDIR}/EP_loops.bedpe ${OUTDIR}/PE_loops.bedpe > ${OUTDIR}/EP_PE_loops.bedpe
cat ${OUTDIR}/EP_PE_loops.bedpe ${OUTDIR}/PP_loops.bedpe > ${OUTDIR}/EP_PE_PP_loops.bedpe
cat ${OUTDIR}/EE_loops.bedpe ${OUTDIR}/EP_PE_PP_loops.bedpe > ${OUTDIR}/all_loops.bedpe

awk -v OFS="\t" '{print $1":"$2":"$3":"$4":"$5":"$6,$1,$2,$3,$4,$5,$6,$7,$8,$9}' ${OUTDIR}/EE_loops.bedpe | sort -u -k1,1 | cut -f 2-8 > ${OUTDIR}/EE_loops_uniq.bedpe 
awk -v OFS="\t" '{print $1":"$2":"$3":"$4":"$5":"$6,$1,$2,$3,$4,$5,$6,$7,$8,$9}' ${OUTDIR}/PP_loops.bedpe | sort -u -k1,1 | cut -f 2-8 > ${OUTDIR}/PP_loops_uniq.bedpe
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
rm -f ${OUTDIR}/a1.bed ${OUTDIR}/a2.bed ${OUTDIR}/*_intersect.txt ${OUTDIR}/*_loops.txt ${OUTDIR}/loops_labeled.bedpe K1 K2 T1 T2

echo "Done."
