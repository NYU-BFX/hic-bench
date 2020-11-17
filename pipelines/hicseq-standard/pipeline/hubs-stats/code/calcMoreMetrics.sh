#!/bin/bash
main_dir=$1
k27ac=$2
tss=$3
atac=$4

echo $main_dir
echo $k27ac
echo $tss
echo $atac

module load gtools
module load bedtools

#k27ac="./params/LOUCY_H3K27ac_peaks_wSignal.bed"
#tss="./params/tss_w_nmre_sorted.bed"
#atac="FALSE"
#atac="./params/LOUCY_H3K27ac_peaks_wSignal.bed"
#main_dir="/gpfs/home/rodrij92/LOUCY-HiChIP-aug25/hic-bench/pipelines/hicseq-standard/pipeline/hubs-stats"

	
## Compute EC, EX BACKGROUND and OBSERVED/EXPECTED metrics ##
echo 'EC.bg' >> EC_bg.txt
echo 'EX.bg' >> EX_bg.txt
echo 'EC.oe' >> EC_oe.txt
echo 'EX.oe' >> EX_oe.txt

awk 'NR>1' EP_master.tsv > temp.txt

awk '{				
if ($14 > 0)			
	print 1/$14		
else				
	print $14		
}' temp.txt >> EC_bg.txt


awk '{				
if ($8 > 0)			
	print 1/$8		
else				
	print $8		
}' temp.txt >> EX_bg.txt


paste EP_master.tsv EC_bg.txt EX_bg.txt > temp2.txt
awk 'NR>1' temp2.txt > temp3.txt

awk '{				
if ($21 > 0)			
	print $4/$21		
else				
	print "NA"		
}' temp3.txt >> EC_oe.txt


awk '{				
if ($22 > 0)			
	print $5/$22		
else				
	print "NA"		
}' temp3.txt >> EX_oe.txt


### PROMOTER ENHANCER ACTIVITY AND PROMOTER ACCESSIBILITY METRICS ###
if [[ $atac != "FALSE" ]]
then
	## intersect ATAC peaks with ENHANCERS and calculate the MEAN/MAX/SUM ENHANCER ACCESSIBILITY
	#4kb
	bedtools window -w 2000 -a ${main_dir}/${k27ac} -b ${main_dir}/${atac} | cut -f5,9 | tools-mergeuniq -merge | tools-vectors m > MEA_in
	bedtools window -v -w 2000 -a ${main_dir}/${k27ac} -b ${main_dir}/${atac} | cut -f5,9 | tools-mergeuniq -merge | tools-vectors m > MEA_out
	cat MEA_in MEA_out | sort -k1,1b > mean.enhancer.accessibility.tsv

	bedtools window -w 2000 -a ${main_dir}/${k27ac} -b ${main_dir}/${atac} | cut -f5,9 | tools-mergeuniq -merge | tools-vectors max > MXEA_in
	bedtools window -v -w 2000 -a ${main_dir}/${k27ac} -b ${main_dir}/${atac} | cut -f5,9 | tools-mergeuniq -merge | tools-vectors max > MXEA_out
	cat MXEA_in MXEA_out | sort -k1,1b > max.enhancer.accessibility.tsv

	bedtools window -w 2000 -a ${main_dir}/${k27ac} -b ${main_dir}/${atac} | cut -f5,9 | tools-mergeuniq -merge | tools-vectors sum > SMEA_in
	bedtools window -v -w 2000 -a ${main_dir}/${k27ac} -b ${main_dir}/${atac} | cut -f5,9 | tools-mergeuniq -merge | tools-vectors sum > SMEA_out
	cat SMEA_in SMEA_out | sort -k1,1b > sum.enhancer.accessibility.tsv

	## intersect ATAC peaks with PROMOTERS and calculate the MEAN/MAX/SUM PROMOTER ACCESSIBILITY
	#4kb
	bedtools window -w 2000 -a ${main_dir}/${tss} -b ${main_dir}/${atac} | cut -f4,10 | tools-mergeuniq -merge | tools-vectors m > MPA_in
	bedtools window -v -w 2000 -a ${main_dir}/${tss} -b ${main_dir}/${atac} | cut -f4,10 | tools-mergeuniq -merge | tools-vectors m > MPA_out
	cat MPA_in MPA_out | sort -k1,1b > mean.promoter.accessibility.tsv

	bedtools window -w 2000 -a ${main_dir}/${tss} -b ${main_dir}/${atac} | cut -f4,10 | tools-mergeuniq -merge | tools-vectors max > MPXA_in
	bedtools window -v -w 2000 -a ${main_dir}/${tss} -b ${main_dir}/${atac} | cut -f4,10 | tools-mergeuniq -merge | tools-vectors max > MPXA_out
	cat MPXA_in MPXA_out | sort -k1,1b > max.promoter.accessibility.tsv

	bedtools window -w 2000 -a ${main_dir}/${tss} -b ${main_dir}/${atac} | cut -f4,10 | tools-mergeuniq -merge | tools-vectors sum > SMPA_in
	bedtools window -v -w 2000 -a ${main_dir}/${tss} -b ${main_dir}/${atac} | cut -f4,10 | tools-mergeuniq -merge | tools-vectors sum > SMPA_out
	cat SMPA_in SMPA_out | sort -k1,1b > sum.promoter.accessibility.tsv
	
	#20kb
	bedtools window -w 10000 -a ${main_dir}/${tss} -b ${main_dir}/${atac} | cut -f4,10 | tools-mergeuniq -merge | tools-vectors m > MPA_in_20
	bedtools window -v -w 10000 -a ${main_dir}/${tss} -b ${main_dir}/${atac} | cut -f4,10 | tools-mergeuniq -merge | tools-vectors m > MPA_out_20
	cat MPA_in_20 MPA_out_20 | sort -k1,1b > mean.promoter.accessibility_20.tsv

	bedtools window -w 10000 -a ${main_dir}/${tss} -b ${main_dir}/${atac} | cut -f4,10 | tools-mergeuniq -merge | tools-vectors max > MPXA_in_20
	bedtools window -v -w 10000 -a ${main_dir}/${tss} -b ${main_dir}/${atac} | cut -f4,10 | tools-mergeuniq -merge | tools-vectors max > MPXA_out_20
	cat MPXA_in_20 MPXA_out_20 | sort -k1,1b > max.promoter.accessibility_20.tsv

	bedtools window -w 10000 -a ${main_dir}/${tss} -b ${main_dir}/${atac} | cut -f4,10 | tools-mergeuniq -merge | tools-vectors sum > SMPA_in_20
	bedtools window -v -w 10000 -a ${main_dir}/${tss} -b ${main_dir}/${atac} | cut -f4,10 | tools-mergeuniq -merge | tools-vectors sum > SMPA_out_20
	cat SMPA_in_20 SMPA_out_20 | sort -k1,1b > sum.promoter.accessibility_20.tsv
	
	#50kb
	bedtools window -w 25000 -a ${main_dir}/${tss} -b ${main_dir}/${atac} | cut -f4,10 | tools-mergeuniq -merge | tools-vectors m > MPA_in_50
	bedtools window -v -w 25000 -a ${main_dir}/${tss} -b ${main_dir}/${atac} | cut -f4,10 | tools-mergeuniq -merge | tools-vectors m > MPA_out_50
	cat MPA_in_50 MPA_out_50 | sort -k1,1b > mean.promoter.accessibility_50.tsv

	bedtools window -w 25000 -a ${main_dir}/${tss} -b ${main_dir}/${atac} | cut -f4,10 | tools-mergeuniq -merge | tools-vectors max > MPXA_in_50
	bedtools window -v -w 25000 -a ${main_dir}/${tss} -b ${main_dir}/${atac} | cut -f4,10 | tools-mergeuniq -merge | tools-vectors max > MPXA_out_50
	cat MPXA_in_50 MPXA_out_50 | sort -k1,1b > max.promoter.accessibility_50.tsv

	bedtools window -w 25000 -a ${main_dir}/${tss} -b ${main_dir}/${atac} | cut -f4,10 | tools-mergeuniq -merge | tools-vectors sum > SMPA_in_50
	bedtools window -v -w 25000 -a ${main_dir}/${tss} -b ${main_dir}/${atac} | cut -f4,10 | tools-mergeuniq -merge | tools-vectors sum > SMPA_out_50
	cat SMPA_in_50 SMPA_out_50 | sort -k1,1b > sum.promoter.accessibility_50.tsv
fi

# intersect K27AC peaks with PROMOTERS and calculate the MEAN/MAX/SUM PROMOTER ENHANCER ACTIVITY
bedtools window -w 2000 -a ${main_dir}/${tss} -b ${main_dir}/${k27ac} | cut -f4,10 | tools-mergeuniq -merge | tools-vectors m  > MPEA_in
bedtools window -v -w 2000 -a ${main_dir}/${tss} -b ${main_dir}/${k27ac} | cut -f4,10 | tools-mergeuniq -merge | tools-vectors m > MPEA_out
cat MPEA_in MPEA_out | sort -k1,1b > mean.promoter.enh.activity.tsv

bedtools window -w 2000 -a ${main_dir}/${tss} -b ${main_dir}/${k27ac} | cut -f4,10 | tools-mergeuniq -merge | tools-vectors max > MXPEA_in
bedtools window -v -w 2000 -a ${main_dir}/${tss} -b ${main_dir}/${k27ac} | cut -f4,10 | tools-mergeuniq -merge | tools-vectors max > MXPEA_out
cat MXPEA_in MXPEA_out | sort -k1,1b > max.promoter.enh.activity.tsv

bedtools window -w 2000 -a ${main_dir}/${tss} -b ${main_dir}/${k27ac} | cut -f4,10 | tools-mergeuniq -merge | tools-vectors sum > SMPEA_in
bedtools window -v -w 2000 -a ${main_dir}/${tss} -b ${main_dir}/${k27ac} | cut -f4,10 | tools-mergeuniq -merge | tools-vectors sum > SMPEA_out
cat SMPEA_in SMPEA_out | sort -k1,1b > sum.promoter.enh.activity.tsv

bedtools window -w 10000 -a ${main_dir}/${tss} -b ${main_dir}/${k27ac} | cut -f4,10 | tools-mergeuniq -merge | tools-vectors m  > MPEA_in_20
bedtools window -v -w 10000 -a ${main_dir}/${tss} -b ${main_dir}/${k27ac} | cut -f4,10 | tools-mergeuniq -merge | tools-vectors m > MPEA_out_20
cat MPEA_in_20 MPEA_out_20 | sort -k1,1b > mean.promoter.enh.activity_20.tsv

bedtools window -w 10000 -a ${main_dir}/${tss} -b ${main_dir}/${k27ac} | cut -f4,10 | tools-mergeuniq -merge | tools-vectors max > MXPEA_in_20
bedtools window -v -w 10000 -a ${main_dir}/${tss} -b ${main_dir}/${k27ac} | cut -f4,10 | tools-mergeuniq -merge | tools-vectors max > MXPEA_out_20
cat MXPEA_in_20 MXPEA_out_20 | sort -k1,1b > max.promoter.enh.activity_20.tsv

bedtools window -w 10000 -a ${main_dir}/${tss} -b ${main_dir}/${k27ac} | cut -f4,10 | tools-mergeuniq -merge | tools-vectors sum > SMPEA_in_20
bedtools window -v -w 10000 -a ${main_dir}/${tss} -b ${main_dir}/${k27ac} | cut -f4,10 | tools-mergeuniq -merge | tools-vectors sum > SMPEA_out_20
cat SMPEA_in_20 SMPEA_out_20 | sort -k1,1b > sum.promoter.enh.activity_20.tsv

bedtools window -w 25000 -a ${main_dir}/${tss} -b ${main_dir}/${k27ac} | cut -f4,10 | tools-mergeuniq -merge | tools-vectors m  > MPEA_in_50
bedtools window -v -w 25000 -a ${main_dir}/${tss} -b ${main_dir}/${k27ac} | cut -f4,10 | tools-mergeuniq -merge | tools-vectors m > MPEA_out_50
cat MPEA_in_50 MPEA_out_50 | sort -k1,1b > mean.promoter.enh.activity_50.tsv

bedtools window -w 25000 -a ${main_dir}/${tss} -b ${main_dir}/${k27ac} | cut -f4,10 | tools-mergeuniq -merge | tools-vectors max > MXPEA_in_50
bedtools window -v -w 25000 -a ${main_dir}/${tss} -b ${main_dir}/${k27ac} | cut -f4,10 | tools-mergeuniq -merge | tools-vectors max > MXPEA_out_50
cat MXPEA_in_50 MXPEA_out_50 | sort -k1,1b > max.promoter.enh.activity_50.tsv

bedtools window -w 25000 -a ${main_dir}/${tss} -b ${main_dir}/${k27ac} | cut -f4,10 | tools-mergeuniq -merge | tools-vectors sum > SMPEA_in_50
bedtools window -v -w 25000 -a ${main_dir}/${tss} -b ${main_dir}/${k27ac} | cut -f4,10 | tools-mergeuniq -merge | tools-vectors sum > SMPEA_out_50
cat SMPEA_in_50 SMPEA_out_50 | sort -k1,1b > sum.promoter.enh.activity_50.tsv


## Add PROMOTER, ENHANCER and ACCESSIBILITY info, and EC, EX OBSERVED/EXPECTED RATIOS ##
sort -k4,4b ${main_dir}/${tss} | cut -f 2-6 > tss_sort.bed
paste temp2.txt EC_oe.txt EX_oe.txt | sort -k2,2b > temp4.txt    # EC EX OE DATA 
join -1 2 -2 3 temp4.txt tss_sort.bed | sort -k3,3b > temp5.txt  # PROMOTER DATA
sort -k5,5b ${main_dir}/${k27ac} | cut -f 2-5 > k27ac_sort.bed

join -1 3 -2 4 temp5.txt k27ac_sort.bed > temp6.txt      # ENHANCER DATA
	
if [[ $atac != "FALSE" ]] # COMPUTE/ADD MEAN ENHANCER ACCESSIBILITY
then		  		 
	join -1 1 -2 1 temp6.txt mean.enhancer.accessibility.tsv | tr ' ' '\t' > temp; mv temp temp6.txt      

else			 # If there is NO ATAC DATA add "1"
	awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$31,$32,$33,$34,$35,1}' temp6.txt > temp; mv temp temp6.txt
fi


## Compute PROMOTER MAX/SUM/MEAN metrics (p/e/pe/pp HUBNESS,ENH.ACTIVITY,EC,EX,LOOP ACTIVITY, IMPACT, p/e/pe INTERACTIVITY ##
awk -v OFS="\t" '{print $2,$1,$4}' temp6.txt | grep '^TSS_' | grep 'ENH_' | cut -f1,3 | sort | tools-mergeuniq -merge | tools-vectors max -n 3 | sort -k1b,1 > p.max.EC.tsv
awk -v OFS="\t" '{print $2,$1,$5}' temp6.txt | grep '^TSS_' | grep 'ENH_' | cut -f1,3 | sort | tools-mergeuniq -merge | tools-vectors max -n 3 | sort -k1b,1 > p.max.EX.tsv
awk -v OFS="\t" '{print $2,$1,$24}' temp6.txt | grep '^TSS_' | grep 'ENH_' | cut -f1,3 | sort | tools-mergeuniq -merge | tools-vectors max -n 3 | sort -k1b,1 > p.max.loop.activity.tsv

awk -v OFS="\t" '{print $2,$1,$35}' temp6.txt | grep '^TSS_' | grep 'ENH_' | cut -f1,3 | sort | tools-mergeuniq -merge | tools-vectors max -n 3 | sort -k1b,1 > p.max.enh.activity.tsv
awk -v OFS="\t" '{print $2,$1,$35}' temp6.txt | grep '^TSS_' | grep 'ENH_' | cut -f1,3 | sort | tools-mergeuniq -merge | tools-vectors sum -n 3 | sort -k1b,1 > p.sum.enh.activity.tsv

awk -v OFS="\t" '{print $2,$1,$8}' temp6.txt | grep '^TSS_' | grep 'ENH_' | cut -f1,3 | sort | tools-mergeuniq -merge | tools-vectors max -n 3 | sort -k1b,1 > p.max.e.hub.tsv
awk -v OFS="\t" '{print $2,$1,$8}' temp6.txt | grep '^TSS_' | grep 'ENH_' | cut -f1,3 | sort | tools-mergeuniq -merge | tools-vectors sum -n 3 | sort -k1b,1 > p.sum.e.hub.tsv

awk -v OFS="\t" '{print $2,$1,$10}' temp6.txt | grep '^TSS_' | grep 'ENH_' | cut -f1,3 | sort | tools-mergeuniq -merge | tools-vectors max -n 3 | sort -k1b,1 > p.max.e.impact.tsv
awk -v OFS="\t" '{print $2,$1,$10}' temp6.txt | grep '^TSS_' | grep 'ENH_' | cut -f1,3 | sort | tools-mergeuniq -merge | tools-vectors sum -n 3 | sort -k1b,1 > p.sum.e.impact.tsv
awk -v OFS="\t" '{print $2,$1,$10}' temp6.txt | grep '^TSS_' | grep 'ENH_' | cut -f1,3 | sort | tools-mergeuniq -merge | tools-vectors m -n 3 | sort -k1b,1 > p.mean.e.impact.tsv

awk -v OFS="\t" '{print $2,$1,$9}' temp6.txt | grep '^TSS_' | grep 'ENH_' | cut -f1,3 | sort | tools-mergeuniq -merge | tools-vectors max -n 3 | sort -k1b,1 > p.max.e.interactivity.tsv
awk -v OFS="\t" '{print $2,$1,$9}' temp6.txt | grep '^TSS_' | grep 'ENH_' | cut -f1,3 | sort | tools-mergeuniq -merge | tools-vectors sum -n 3 | sort -k1b,1 > p.sum.e.interactivity.tsv

awk -v OFS="\t" '{print $2,$1,$13}' temp6.txt | grep '^TSS_' | grep 'ENH_' | cut -f1,3 | sort | tools-mergeuniq -merge | tools-vectors max -n 3 | sort -k1b,1 > p.max.e.p.interactivity.tsv
awk -v OFS="\t" '{print $2,$1,$13}' temp6.txt | grep '^TSS_' | grep 'ENH_' | cut -f1,3 | sort | tools-mergeuniq -merge | tools-vectors sum -n 3 | sort -k1b,1 > p.sum.e.p.interactivity.tsv

awk -v OFS="\t" '{print $2,$1,$17}' temp6.txt | grep '^TSS_' | grep 'ENH_' | cut -f1,3 | sort | tools-mergeuniq -merge | tools-vectors max -n 3 | sort -k1b,1 > p.max.p.p.interactivity.tsv
awk -v OFS="\t" '{print $2,$1,$17}' temp6.txt | grep '^TSS_' | grep 'ENH_' | cut -f1,3 | sort | tools-mergeuniq -merge | tools-vectors sum -n 3 | sort -k1b,1 > p.sum.p.p.interactivity.tsv

awk -v OFS="\t" '{print $2,$1,$12}' temp6.txt | grep '^TSS_' | grep 'ENH_' | cut -f1,3 | sort | tools-mergeuniq -merge | tools-vectors max -n 3 | sort -k1b,1 > p.max.e.p.hub.tsv
awk -v OFS="\t" '{print $2,$1,$12}' temp6.txt | grep '^TSS_' | grep 'ENH_' | cut -f1,3 | sort | tools-mergeuniq -merge | tools-vectors sum -n 3 | sort -k1b,1 > p.sum.e.p.hub.tsv

awk -v OFS="\t" '{print $2,$1,$14}' temp6.txt | grep '^TSS_' | grep 'ENH_' | cut -f1,3 | sort | tools-mergeuniq -merge | tools-vectors max -n 3 | sort -k1b,1 > max.p.hub.tsv
awk -v OFS="\t" '{print $2,$1,$14}' temp6.txt | grep '^TSS_' | grep 'ENH_' | cut -f1,3 | sort | tools-mergeuniq -merge | tools-vectors sum -n 3 | sort -k1b,1 > sum.p.hub.tsv

awk -v OFS="\t" '{print $2,$1,$15}' temp6.txt | grep '^TSS_' | grep 'ENH_' | cut -f1,3 | sort | tools-mergeuniq -merge | tools-vectors max -n 3 | sort -k1b,1 > max.p.p.hub.tsv
awk -v OFS="\t" '{print $2,$1,$15}' temp6.txt | grep '^TSS_' | grep 'ENH_' | cut -f1,3 | sort | tools-mergeuniq -merge | tools-vectors sum -n 3 | sort -k1b,1 > sum.p.p.hub.tsv

awk -v OFS="\t" '{print $2,$1,$16}' temp6.txt | grep '^TSS_' | grep 'ENH_' | cut -f1,3 | sort | tools-mergeuniq -merge | tools-vectors max -n 3 | sort -k1b,1 > max.p.e.hub.tsv
awk -v OFS="\t" '{print $2,$1,$16}' temp6.txt | grep '^TSS_' | grep 'ENH_' | cut -f1,3 | sort | tools-mergeuniq -merge | tools-vectors sum -n 3 | sort -k1b,1 > sum.p.e.hub.tsv


### 36 columns until here (temp6) ###

## Compute abc-metrics ##
awk '{print ($35+$36)/2}' temp6.txt > mean.enh.atac_activity.txt                # mean.enh.atac_activity
paste temp6.txt mean.enh.atac_activity.txt > t; mv t temp6.txt

awk -v OFS="\t" '{print $1,$2,$3,$37*$24}' temp6.txt > abc.tsv                  # abc value
cut -f2,3,4 abc.tsv | grep '^TSS_' | grep 'ENH_' | cut -f1,3 | sort | tools-mergeuniq -merge | tools-vectors sum -n 3 | sort -k1b,1 > sum.p.abc.tsv         # sum.p.abc value 

# compute enhancer centric ABC metrics
cut -f2,4 abc.tsv | sort | tools-mergeuniq -merge | tools-vectors sum -n 3 > e.sum.abc.tsv
cut -f2,4 abc.tsv | sort | tools-mergeuniq -merge | tools-vectors m -n 3 > e.mean.abc.tsv
cut -f2,4 abc.tsv | sort | tools-mergeuniq -merge | tools-vectors max -n 3 > e.max.abc.tsv
cut -f 3,4 abc.tsv | sort -k1b,1 > t; mv t abc.tsv                                   # abc value clean
join e.sum.abc.tsv e.mean.abc.tsv | join - e.max.abc.tsv | tr ' ' '\t' | sort -r -n -k2,2 >> e.abc_metrics.tsv

# gene annotation of enhancer centric ABC metrics
echo "enhancer.id e.sum.abc e.mean.abc e.max.abc gene.id EP.id" | tr ' ' '\t' >> e.abc_metrics_geneAnnot.tsv
cut -f 1-3 temp6.txt | sort -k3,3b | tr '\t' ' ' |  tools-cols 1 0 2 | sort -k1,1b > EP.tsv
awk 'NR>1' e.abc_metrics.tsv | sort -k1,1b | join - EP.tsv | sort -r -n -k2,2 | tr ' ' '\t' >> e.abc_metrics_geneAnnot.tsv	

### 37 columns until here (temp6) ###

#### MAKE FINAL MATRIX ####
sort -k3,3b temp6.txt > temp; mv temp temp6.txt

if [[ $atac != "FALSE" ]]	# WITH ATAC DATA
then
	join -1 3 -2 1 temp6.txt abc.tsv | sort -k2,2b | join -1 2 -2 1 - max.enhancer.accessibility.tsv  | join - sum.enhancer.accessibility.tsv | tr ' ' '\t' | sort -k3,3b > temp; mv temp temp6.txt
	### 40 columns until here (temp6) ###

	join -1 3 -2 1 temp6.txt sum.p.abc.tsv  | join - p.max.EC.tsv | join - p.max.EX.tsv | join - p.max.loop.activity.tsv | join - p.max.enh.activity.tsv | join - p.sum.enh.activity.tsv | join - p.max.e.hub.tsv | join - p.sum.e.hub.tsv | join - p.max.e.impact.tsv | join - p.sum.e.impact.tsv | join - p.mean.e.impact.tsv | join - p.max.e.interactivity.tsv | join - p.sum.e.interactivity.tsv | join - p.max.e.p.interactivity.tsv | join - p.sum.e.p.interactivity.tsv | join - p.max.p.p.interactivity.tsv | join - p.sum.p.p.interactivity.tsv | join - p.max.e.p.hub.tsv | join - p.sum.e.p.hub.tsv | join - max.p.hub.tsv | join - sum.p.hub.tsv | join - max.p.p.hub.tsv | join - sum.p.p.hub.tsv | join - max.p.e.hub.tsv | join - sum.p.e.hub.tsv | tr ' ' '\t' > t; mv t temp6.txt
	
	# add abc_score
	awk '{				
	if ($41 > 0)			
		print $38/$41		
	else				
		print $41		
	}' temp6.txt > abc_score.txt

	paste temp6.txt abc_score.txt > t; mv t temp6.txt
	
	### 66 columns until here (temp6) ###

	join temp6.txt mean.promoter.enh.activity.tsv | join - max.promoter.enh.activity.tsv | join - sum.promoter.enh.activity.tsv | join - mean.promoter.enh.activity_20.tsv | join - max.promoter.enh.activity_20.tsv | join - sum.promoter.enh.activity_20.tsv | join - mean.promoter.enh.activity_50.tsv | join - max.promoter.enh.activity_50.tsv | join - sum.promoter.enh.activity_50.tsv | join - mean.promoter.accessibility.tsv | join - max.promoter.accessibility.tsv | join - sum.promoter.accessibility.tsv | join - mean.promoter.accessibility_20.tsv | join - max.promoter.accessibility_20.tsv | join - sum.promoter.accessibility_20.tsv | join - mean.promoter.accessibility_50.tsv | join - max.promoter.accessibility_50.tsv | join - sum.promoter.accessibility_50.tsv | tr ' ' '\t' > t; mv t temp6.txt
	### 84 columns until here (temp6) ###

	# FINAL MATRIX (WITH ATAC)
	echo "EP.id promoter.id enhancer.id EC EX EC.EX_sum EC.EX_diff e.hub e.interactivity e.impact p.EX.sum e.phub e.p.interactivity p.hub pp.hub pe.hub pp.interactivity pe.interactivity p.interactivity loop.chr loop.start loop.end loop.distance loop.activity EC.bg EX.bg EC.oe EX.oe tss_start tss_end tss_ensembleID tss_strand enh_start enh_end enh_signal mean.enh.accessibility mean.enh.atac.activity abc max.enh.accessibility sum.enh.accessibility sum.p.abc p.max.EC p.max.EX p.max.loop.activity p.max.enh.activity p.sum.enh.activity p.max.e.hub p.sum.e.hub p.max.e.impact p.sum.e.impact p.mean.e.impact p.max.e.interactivity p.sum.e.interactivity p.max.e.p.interactivity p.sum.e.p.interactivity p.max.p.p.interactivity p.sum.p.p.interactivity p.max.e.p.hub p.sum.e.p.hub max.p.hub sum.p.hub max.p.p.hub sum.p.p.hub max.p.e.hub sum.p.e.hub abc_score mean.promoter.enh.activity max.promoter.enh.activity sum.promoter.enh.activity mean.promoter.enh.activity_20 max.promoter.enh.activity_20 sum.promoter.enh.activity_20 mean.promoter.enh.activity_50 max.promoter.enh.activity_50 sum.promoter.enh.activity_50 mean.promoter.accessibility max.promoter.accessibility sum.promoter.accessibility mean.promoter.accessibility_20 max.promoter.accessibility_20 sum.promoter.accessibility_20 mean.promoter.accessibility_50 max.promoter.accessibility_50 sum.promoter.accessibility_50" | tr ' ' '\t' >> EP_metrics_full.tsv 
	
	awk -v OFS="\t" '{print $3,$1,$2}' temp6.txt > chunk1.tsv
	cut -f4-84 temp6.txt > chunk2.tsv
	paste chunk1.tsv chunk2.tsv >> EP_metrics_full.tsv

	# make promoter centric table
	sort -u -k2,2b  EP_metrics_full.tsv | cut -f2,11,14-19,41-84 > p.centric_metrics.tsv


else	### IF NO ATAC DATA ###
	join -1 3 -2 1 temp6.txt abc.tsv | tr ' ' '\t' | sort -k3,3b > temp; mv temp temp6.txt
	### 38 columns until here (temp6) ###

	join -1 3 -2 1 temp6.txt sum.p.abc.tsv  | join - p.max.EC.tsv | join - p.max.EX.tsv | join - p.max.loop.activity.tsv | join - p.max.enh.activity.tsv | join - p.sum.enh.activity.tsv | join - p.max.e.hub.tsv | join - p.sum.e.hub.tsv | join - p.max.e.impact.tsv | join - p.sum.e.impact.tsv | join - p.mean.e.impact.tsv | join - p.max.e.interactivity.tsv | join - p.sum.e.interactivity.tsv | join - p.max.e.p.interactivity.tsv | join - p.sum.e.p.interactivity.tsv | join - p.max.p.p.interactivity.tsv | join - p.sum.p.p.interactivity.tsv | join - p.max.e.p.hub.tsv | join - p.sum.e.p.hub.tsv | join - max.p.hub.tsv | join - sum.p.hub.tsv | join - max.p.p.hub.tsv | join - sum.p.p.hub.tsv | join - max.p.e.hub.tsv | join - sum.p.e.hub.tsv | tr ' ' '\t' > t; mv t temp6.txt
        
	# add abc_score
        awk '{
	if ($39 > 0)
                print $38/$39
        else
            	print $39
        }' temp6.txt > abc_score.txt
	
        paste temp6.txt abc_score.txt > t; mv t temp6.txt

	### 64 columns until here (temp6) ###

	join temp6.txt mean.promoter.enh.activity.tsv | join - max.promoter.enh.activity.tsv | join - sum.promoter.enh.activity.tsv | join - mean.promoter.enh.activity_20.tsv | join - max.promoter.enh.activity_20.tsv | join - sum.promoter.enh.activity_20.tsv | join - mean.promoter.enh.activity_50.tsv | join - max.promoter.enh.activity_50.tsv | join - sum.promoter.enh.activity_50.tsv | tr ' ' '\t' > t; mv t temp6.txt
	### 73 columns until here (temp6) ###

	# FINAL MATRIX (NO ATAC)
	echo "EP.id promoter.id enhancer.id EC EX EC.EX_sum EC.EX_diff e.hub e.interactivity e.impact p.EX.sum e.phub e.p.interactivity p.hub pp.hub pe.hub pp.interactivity pe.interactivity p.interactivity loop.chr loop.start loop.end loop.distance loop.activity EC.bg EX.bg EC.oe EX.oe tss_start tss_end tss_ensembleID tss_strand enh_start enh_end enh_signal mean.enh.accessibility mean.enh.atac.activity abc sum.p.abc p.max.EC p.max.EX p.max.loop.activity p.max.enh.activity p.sum.enh.activity p.max.e.hub p.sum.e.hub p.max.e.impact p.sum.e.impact p.mean.e.impact p.max.e.interactivity p.sum.e.interactivity p.max.e.p.interactivity p.sum.e.p.interactivity p.max.p.p.interactivity p.sum.p.p.interactivity p.max.e.p.hub p.sum.e.p.hub max.p.hub sum.p.hub max.p.p.hub sum.p.p.hub max.p.e.hub sum.p.e.hub abc_score mean.promoter.enh.activity max.promoter.enh.activity sum.promoter.enh.activity mean.promoter.enh.activity_20 max.promoter.enh.activity_20 sum.promoter.enh.activity_20 mean.promoter.enh.activity_50 max.promoter.enh.activity_50 sum.promoter.enh.activity_50" | tr ' ' '\t' >> EP_metrics_full.tsv 

	awk -v OFS="\t" '{print $3,$1,$2}' temp6.txt > chunk1.tsv
	cut -f4-73 temp6.txt > chunk2.tsv
	paste chunk1.tsv chunk2.tsv >> EP_metrics_full.tsv

	# make promoter centric table
	sort -u -k2,2b  EP_metrics_full.tsv | cut -f2,11,14-19,39-73 > p.centric_metrics.tsv

fi

mkdir data
mv *.tsv data
mv abc_score.txt data/
mv data/EP_metrics_full.tsv ./
mv data/p.centric_metrics.tsv ./
mv data/e.abc_metrics.tsv ./
mv data/e.abc_metrics_geneAnnot.tsv ./

mkdir loops
mv bedpe2V5C/*uniq.bedpe loops/

## CLEAN UP ##
rm -f EC_bg.txt EX_bg.txt EC_oe.txt EX_oe.txt temp*.txt tss_sort.bed k27ac_sort.bed chunk* *_in *_out *_in_* *_out_* EP_master_uniq.tsv EP_master.tsv mean.enh.atac_activity.txt *EP.tsv
rm -fr bedpe2V5C
