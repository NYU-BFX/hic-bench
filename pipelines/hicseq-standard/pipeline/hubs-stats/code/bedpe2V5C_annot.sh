#!/bin/bash
BEDPE=$1
K27AC=$2
TSS=$3
ATAC=$4
accessible_only=$5
tss_extension=$6
promoter_k27ac_only=$7
standarize_cpm=$8
use_topLoops=$9
winsize=${10}
OUTDIR=${11}

module load bedtools/2.27.1

rm -fr ${OUTDIR}
mkdir -p ${OUTDIR}

# add loop identifier column
awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$1":"$2":"$3":"$4":"$5":"$6}' ${BEDPE} > ${OUTDIR}/temp; mv ${OUTDIR}/temp ${OUTDIR}/loops_labeled.bedpe

# separate loop anchors
cut -f 1-3,7-8 ${OUTDIR}/loops_labeled.bedpe | awk -v OFS="\t" '{print $1,$2,$3,$1":"$2":"$3,$4,$5}' > ${OUTDIR}/a1.bed
cut -f 4-6,7-8 ${OUTDIR}/loops_labeled.bedpe | awk -v OFS="\t" '{print $1,$2,$3,$1":"$2":"$3,$4,$5}' > ${OUTDIR}/a2.bed

# clean k27ac file
K27AC_original=${K27AC}
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

if [[ $k27ac_in_TSS_anchor = "FALSE" ]]
then
	# remove k27ac peaks that are in a TSS anchor
	cat ${OUTDIR}/k27ac_a1_intersect.txt | sort -k8b,8 > ${OUTDIR}/K1
	cat ${OUTDIR}/tss_a1_intersect.txt | sort -k10b,10 > ${OUTDIR}/T1

	cat ${OUTDIR}/k27ac_a2_intersect.txt | sort -k8b,8 > ${OUTDIR}/K2
	cat ${OUTDIR}/tss_a2_intersect.txt | sort -k10b,10 > ${OUTDIR}/T2

	join -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11  -v 1 -1 8 -2 10 ${OUTDIR}/K1 ${OUTDIR}/T1 | tr ' ' '\t' | sort -k10b,10 > ${OUTDIR}/temp; mv ${OUTDIR}/temp ${OUTDIR}/k27ac_a1_intersect.txt
	join -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11  -v 1 -1 8 -2 10 ${OUTDIR}/K2 ${OUTDIR}/T2 | tr ' ' '\t' | sort -k10b,10 > ${OUTDIR}/temp; mv ${OUTDIR}/temp ${OUTDIR}/k27ac_a2_intersect.txt

	wc -l ${OUTDIR}/k27ac_a1_intersect.txt
	wc -l ${OUTDIR}/k27ac_a2_intersect.txt
fi

# get EE,EP,PE,PP loops
join -j10 ${OUTDIR}/k27ac_a1_intersect.txt ${OUTDIR}/k27ac_a2_intersect.txt | awk '$9!=$19'> ${OUTDIR}/EE_loops.txt
join -j12 ${OUTDIR}/tss_a1_intersect.txt ${OUTDIR}/tss_a2_intersect.txt | awk '$11!=$23' > ${OUTDIR}/PP_loops.txt
join -1 10 -2 12 ${OUTDIR}/k27ac_a1_intersect.txt ${OUTDIR}/tss_a2_intersect.txt | awk '$9!=$21' > ${OUTDIR}/EP_loops.txt
join -1 12 -2 10 ${OUTDIR}/tss_a1_intersect.txt ${OUTDIR}/k27ac_a2_intersect.txt | awk '$11!=$21' > ${OUTDIR}/PE_loops.txt

# get loops with at least one k27ac
bedtools intersect -a ${K27AC} -b ${OUTDIR}/a1.bed -wo | sort -k10b,10 > ${OUTDIR}/k27ac_a1_intersect_all.txt
bedtools intersect -a ${K27AC} -b ${OUTDIR}/a2.bed -wo | sort -k10b,10 > ${OUTDIR}/k27ac_a2_intersect_all.txt

cat ${OUTDIR}/k27ac_a1_intersect_all.txt  ${OUTDIR}/k27ac_a2_intersect_all.txt | sort -u -k10,10b | cut -f 10 > ${OUTDIR}/loops_with_k27ac_ids.txt
sort -u -k8,8b ${OUTDIR}/loops_labeled.bedpe > ${OUTDIR}/t; mv ${OUTDIR}/t ${OUTDIR}/loops_labeled.bedpe
sort -k1,1b ${OUTDIR}/loops_with_k27ac_ids.txt > ${OUTDIR}/t; mv ${OUTDIR}/t ${OUTDIR}/loops_with_k27ac_ids.txt
join -1 8 ${OUTDIR}/loops_labeled.bedpe ${OUTDIR}/loops_with_k27ac_ids.txt | awk -v OFS="\t" '{print $2,$3,$4,$5,$6,$7,$8,$1}' | sort -u -k8,8b | cut -f 1-7 > ${OUTDIR}/loops_with_k27ac_uniq.bedpe

# get loops with k27ac in both anchors (no filters)
join -j10 ${OUTDIR}/k27ac_a1_intersect_all.txt ${OUTDIR}/k27ac_a2_intersect_all.txt | awk '$9!=$19'> ${OUTDIR}/EE_loops_noFilters.txt
awk -v OFS="\t" '{print $6,$7,$8,$16,$17,$18,$10,$5,$15}' ${OUTDIR}/EE_loops_noFilters.txt > ${OUTDIR}/EE_loops_noFilters.bedpe
awk -v OFS="\t" '{print $1":"$2":"$3":"$4":"$5":"$6,$1,$2,$3,$4,$5,$6,$7,$8,$9}' ${OUTDIR}/EE_loops_noFilters.bedpe | sort -u -k1,1 | cut -f 2-8 > ${OUTDIR}/EE_loops_noFilters_uniq.bedpe 

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
cat ${OUTDIR}/all_loops_uniq.bedpe | wc -l

echo "Total number of unique loops with PP, EP or PE annotation:"
cat ${OUTDIR}/EP_PE_PP_loops_uniq.bedpe | wc -l

# get loops in virtual5C.csv format
awk -v OFS="\t" '{print $8":"$9,$8,$9,$1,($2+$3)/2,($5+$6)/2,$5-$2,"-",$7,"-","-","0","0","-","-"}' ${OUTDIR}/all_loops.bedpe | sort -u -k1,1 > ${OUTDIR}/all_loops_v5cFormat.tsv

# add reversed duplicated loops (for v5cFormat consistency)
awk -v OFS="\t" '{print $3":"$2,$3,$2,$4,$6,$5,$7,$8,$9,$10,$11,$12,$13,$14,$15}' ${OUTDIR}/all_loops_v5cFormat.tsv | sort -u -k1,1 > ${OUTDIR}/all_loops_rev_v5cFormat.tsv

# get final v5cFormat loops file
echo "Source.anchor.label,Target.anchor.label,Chromosome,Source.anchor.position,Target.anchor.position,Anchor.distance,Count,CPK2B,VP.total.count,Count.sum,pvalue,fdr,p.expected,p.observed" >> ${OUTDIR}/all_loops_wRev_v5cFormat.csv
cat ${OUTDIR}/all_loops_v5cFormat.tsv ${OUTDIR}/all_loops_rev_v5cFormat.tsv | sort -u -k1,1 | cut -f 2-15 | tr '\t' ',' >> ${OUTDIR}/all_loops_wRev_v5cFormat.csv

head ${OUTDIR}/all_loops_wRev_v5cFormat.csv 
# select top loops
if [[ $use_topLoops != "FALSE" ]]
then
	top=$(( 2*$use_topLoops ))	# this is because for each top loop we have 2 loops in the v5cFoormat file (1 original and 1 reversed)
	echo 'We will only use the top '$use_topLoops' loops'
	echo "Source.anchor.label,Target.anchor.label,Chromosome,Source.anchor.position,Target.anchor.position,Anchor.distance,Count,CPK2B,VP.total.count,Count.sum,pvalue,fdr,p.expected,p.observed" >> ${OUTDIR}/all_loops_wRev_v5cFormat_top.csv
	cat ${OUTDIR}/all_loops_wRev_v5cFormat.csv | tr ',' '\t' | cut -f 3-5,8 | awk 'NR>1' | awk -v OFS="\t" '{print $1,$2,$3,$4,$1":"$2":"$3}' | sort -u -k5b,5 | sort -r -k4,4 | head -n $top | cut -f 5 > ${OUTDIR}/topLoops_ids.txt
	cat ${OUTDIR}/all_loops_wRev_v5cFormat.csv | tr ',' '\t' | awk 'NR>1' | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$3":"$4":"$5}' | sort -r -k15b,15 | grep -F -w -f ${OUTDIR}/topLoops_ids.txt | cut -f 1-14 | tr '\t' ',' >> ${OUTDIR}/all_loops_wRev_v5cFormat_top.csv
	mv ${OUTDIR}/all_loops_wRev_v5cFormat_top.csv ${OUTDIR}/all_loops_wRev_v5cFormat.csv
fi

echo "Total loops obtained (with Rev and dups):"
nld=`wc -l ${OUTDIR}/all_loops_wRev_v5cFormat.csv | awk '{print $1}'`
echo $nld
echo "Total unique loops obtained (with Rev):"
nlu=`cat ${OUTDIR}/all_loops_wRev_v5cFormat.csv | tr ',' '\t' | cut -f 3-5,8 | awk 'NR>1' | awk -v OFS="\t" '{print $1,$2,$3,$4,$1":"$2":"$3}' | sort -u -k5b,5 | wc -l`    
echo $nlu

# standarize cpm values so the mean cpm value is 1
if [[ $standarize_cpm = "TRUE" ]]
then
	mean_cpm=`awk 'NR>1' ${OUTDIR}/all_loops_wRev_v5cFormat.csv | tr ',' '\t'|  awk '$5>$4' | awk -v OFS="\t" '{print $3":"$4":"$5,$8}' | sort -u -k1b,1 | cut -f2 | awk '{ total += $1 } END { print total/NR }'`
	mfc=`echo $mean_cpm | awk '{print 1/$1}'`
	echo 'The original mean cpm value is = '$mean_cpm

	echo 'The correction factor is = '$mfc

	head -n1 ${OUTDIR}/all_loops_wRev_v5cFormat.csv >> ${OUTDIR}/temp.txt

	awk 'NR>1' ${OUTDIR}/all_loops_wRev_v5cFormat.csv | tr ',' '\t' | awk -v OFS="\t" -v f="$mfc" '{print $1,$2,$3,$4,$5,$6,$7,$8,f,$9,$10,$11,$12,$13,$14}' | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$8*$9,$10,$11,$12,$13,$14,$15}' | tr '\t' ',' >> ${OUTDIR}/temp.txt
	mv ${OUTDIR}/temp.txt ${OUTDIR}/all_loops_wRev_v5cFormat.csv

	mean_scaled=`awk 'NR>1' ${OUTDIR}/all_loops_wRev_v5cFormat.csv | tr ',' '\t' | awk '$5>$4' | awk -v OFS="\t" '{print $3":"$4":"$5,$8}' | sort -u -k1b,1 | cut -f2 | awk '{ total += $1 } END { print total/NR }'`	
	echo 'The mean cpm value after correction is = '$mean_scaled
fi

# make final loops bedpe
awk -F"," -v OFS="\t" '{print $3,$4,$5,$8}' ${OUTDIR}/all_loops_wRev_v5cFormat.csv | awk 'NR>1' > ${OUTDIR}/topLoops.tsv

awk ' BEGIN { OFS="\t"} {                             \
	if ($2 < $3)                                  \
                print $1,$2,$3,$4;                       \
        else                                          \
                print $1,$3,$2,$4;                       \
}'  ${OUTDIR}/topLoops.tsv | awk -v OFS="\t" '{print $1,$2,$3,$4,$1":"$2":"$3}' | sort -u -k5,5b | cut -f 1-4 > ${OUTDIR}/topLoops.bed

awk -v OFS="\t" -v ext=$winsize '{print $1,$2-(ext/2),$2+(ext/2),$1,$3-(ext/2),$3+(ext/2),$4}' ${OUTDIR}/topLoops.bed > ${OUTDIR}/topLoops.bedpe
rm -f ${OUTDIR}/topLoops.tsv ${OUTDIR}/topLoops.bed

### ADD XP, PX and XX loop data ###

# anchors with no k27ac
bedtools intersect -v -a ${OUTDIR}/a1.bed -b ${OUTDIR}/k27ac_clean.bed -wo | sort -k6b,6 > ${OUTDIR}/a1_no_k27ac.txt
bedtools intersect -v -a ${OUTDIR}/a2.bed -b ${OUTDIR}/k27ac_clean.bed -wo | sort -k6b,6 > ${OUTDIR}/a2_no_k27ac.txt

# anchors with no k27ac and no tss
bedtools intersect -v -a ${OUTDIR}/a1_no_k27ac.txt -b ${TSS} -wo | sort -k6b,6 > ${OUTDIR}/a1_no_k27ac_tss.txt
bedtools intersect -v -a ${OUTDIR}/a2_no_k27ac.txt -b ${TSS} -wo | sort -k6b,6 > ${OUTDIR}/a2_no_k27ac_tss.txt

if [[ $ATAC != "FALSE" ]]
then
	# anchors without k27ac and tss, with ATAC ("X")
	bedtools intersect -a ${OUTDIR}/a1_no_k27ac_tss.txt -b ${ATAC} -wo | sort -k6b,6 | cut -f 1-6 > ${OUTDIR}/a1_X.txt
	bedtools intersect -a ${OUTDIR}/a2_no_k27ac_tss.txt -b ${ATAC} -wo | sort -k6b,6 | cut -f 1-6 > ${OUTDIR}/a2_X.txt
else
	cp ${OUTDIR}/a1_no_k27ac_tss.txt ${OUTDIR}/a1_X.txt
        cp ${OUTDIR}/a2_no_k27ac_tss.txt ${OUTDIR}/a2_X.txt
fi

# loops whose both anchors are "X" (-> XX_loops.txt)
join -1 6 -2 6 ${OUTDIR}/a1_X.txt ${OUTDIR}/a2_X.txt | awk -v OFS="\t" '{print $2,$3,$4,$7,$8,$9,$6,$1,$5,$10}' | sort -u -k8,8b > ${OUTDIR}/XX_loops.txt

# get XX bedpe and XX loop ids
cut -f 1-7 ${OUTDIR}/XX_loops.txt > ${OUTDIR}/XX_loops_uniq.bedpe
cut -f8 ${OUTDIR}/XX_loops.txt > ${OUTDIR}/XX_loop_ids.txt

# get XP_PX loops (.txt)
cut -f8 ${OUTDIR}/loops_labeled.bedpe > ${OUTDIR}/all_loop_ids.txt
grep -v -x -F -f ${OUTDIR}/XX_loop_ids.txt ${OUTDIR}/all_loop_ids.txt | sort -k1,1b > ${OUTDIR}/XP_PX_loop_ids.txt

join -1 8 ${OUTDIR}/loops_labeled.bedpe ${OUTDIR}/XP_PX_loop_ids.txt | awk -v OFS="\t" '{print $2,$3,$4,$5,$6,$7,$8,$1}' | sort -u -k8,8b > ${OUTDIR}/XP_PX_loops.txt

# get PX and XP loops separately
# get PX and XP loops
bedtools intersect -a ${OUTDIR}/a1.bed -b ${TSS} -wa | cut -f 1-3,6 | sort -k4,4b > ${OUTDIR}/a1_tss.txt
bedtools intersect -a ${OUTDIR}/a2.bed -b ${TSS} -wa | cut -f 1-3,6 | sort -k4,4b > ${OUTDIR}/a2_tss.txt
bedtools intersect -v -a ${OUTDIR}/a1.bed -b ${TSS} -wa | cut -f 1-3,6 | bedtools intersect -v -a stdin -b ${K27AC} | sort -k4,4b > ${OUTDIR}/a1_no_tss_enh.txt
bedtools intersect -v -a ${OUTDIR}/a2.bed -b ${TSS} -wa | cut -f 1-3,6 | bedtools intersect -v -a stdin -b ${K27AC} | sort -k4,4b > ${OUTDIR}/a2_no_tss_enh.txt

if [[ $ATAC != "FALSE" ]]
then
	bedtools intersect -a ${OUTDIR}/a1_no_tss_enh -b ${ATAC} | sort -k4,4b > ${OUTDIR}/a1_X.txt
	bedtools intersect -a ${OUTDIR}/a2_no_tss_enh -b ${ATAC} | sort -k4,4b > ${OUTDIR}/a2_X.txt
else
	cp ${OUTDIR}/a1_no_tss_enh.txt ${OUTDIR}/a1_X.txt 	
	cp ${OUTDIR}/a2_no_tss_enh.txt ${OUTDIR}/a2_X.txt 	
fi

join -1 4 -2 4 ${OUTDIR}/a1_tss.txt ${OUTDIR}/a2_X.txt | sort -u -k1,1b | cut -f1 > ${OUTDIR}/PX_ids.txt
join -1 4 -2 4 ${OUTDIR}/a1_X.txt ${OUTDIR}/a2_tss.txt | sort -u -k1,1b | cut -f1 > ${OUTDIR}/XP_ids.txt

cat ${OUTDIR}/PX_ids.txt ${OUTDIR}/XP_ids.txt | sort -u -k1,1b | cut -d " " -f1 | sort > ${OUTDIR}/PX_XP_ids.txt

sort -k8,8b ${OUTDIR}/loops_labeled.bedpe > ${OUTDIR}/temp; mv ${OUTDIR}/temp ${OUTDIR}/loops_labeled.bedpe
join -1 8 -2 1  ${OUTDIR}/loops_labeled.bedpe ${OUTDIR}/PX_XP_ids.txt | tr ' ' '\t' | cut -f 2-8 > ${OUTDIR}/PX_XP_loops_uniq.bedpe
cut -f 1-3,7 ${OUTDIR}/PX_XP_loops_uniq.bedpe > ${OUTDIR}/PX_XP_a1.txt
cut -f 4-7 ${OUTDIR}/PX_XP_loops_uniq.bedpe > ${OUTDIR}/PX_XP_a2.txt

bedtools intersect -a ${TSS} -b ${OUTDIR}/PX_XP_a1.txt -wa -wb | cut -f 4,10 > ${OUTDIR}/PX.tsv
bedtools intersect -a ${TSS} -b ${OUTDIR}/PX_XP_a2.txt -wa -wb | cut -f 4,10 > ${OUTDIR}/XP.tsv
cat ${OUTDIR}/PX.tsv ${OUTDIR}/XP.tsv > ${OUTDIR}/PX_XP.tsv

cut -f1 ${OUTDIR}/PX_XP.tsv | sort | uniq -c | tools-cols 1 0 | tr ' ' '\t' > ${OUTDIR}/p.x.hubness.tsv
cat ${OUTDIR}/PX_XP.tsv | sort | tools-mergeuniq -merge | tools-vectors sum -n 3 > ${OUTDIR}/p.x.sum.interactivity.tsv
cat ${OUTDIR}/PX_XP.tsv | sort | tools-mergeuniq -merge | tools-vectors m -n 3 > ${OUTDIR}/p.x.mean.interactivity.tsv
cat ${OUTDIR}/PX_XP.tsv | sort | tools-mergeuniq -merge | tools-vectors max -n 3 > ${OUTDIR}/p.x.max.interactivity.tsv
echo "gene.id p.x.hubness p.x.sum.interactivity p.x.mean.interactivity p.x.max.interactivity" | tr ' ' '\t' >> ${OUTDIR}/p.x.metrics.tsv
join  ${OUTDIR}/p.x.hubness.tsv ${OUTDIR}/p.x.sum.interactivity.tsv | join - ${OUTDIR}/p.x.mean.interactivity.tsv | join - ${OUTDIR}/p.x.max.interactivity.tsv | tr ' ' '\t' >> ${OUTDIR}/p.x.metrics.tsv 

### ADD k27ac SIGNAL TO LOOPS ###

./code/scripts-hubs-stats-addSignal2Loops.sh loops_with_k27ac_uniq.bedpe ${K27AC_original} ${OUTDIR} loops_with_k27ac_uniq.tsv 
./code/scripts-hubs-stats-addSignal2Loops.sh EE_loops_noFilters_uniq.bedpe ${K27AC_original} ${OUTDIR} EE_loops_noFilters_uniq.tsv
./code/scripts-hubs-stats-addSignal2Loops.sh loops_labeled.bedpe ${K27AC_original} ${OUTDIR} loops_labeled_k27ac_signal.tsv

# clean up
rm -f ${OUTDIR}/a1.bed ${OUTDIR}/a2.bed ${OUTDIR}/*_intersect.txt ${OUTDIR}/*_loops.txt K1 K2 T1 T2
rm -f ${OUTDIR}/a1_X* ${OUTDIR}/a2_X* ${OUTDIR}/XP_PX_loops.txt ${OUTDIR}/*ids.txt ${OUTDIR}/a1_no_* ${OUTDIR}/a2_no_*

echo "Done."
