#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: hicseq-template.tcsh OUTPUT-DIR PARAM-SCRIPT BRANCH OBJECT(S)
##

if ($#argv != 4) then
  grep '^##' $0
  exit
endif

set outdir = $1
set params = $2
set branch = $3
set objects = ($4)

# does this operation allow grouping of input objects?
if ($#objects>1) then
  send2err "Error: this operation is not implemented for multi-object grouping."
  exit 1
else 
  set object = $objects[1]
endif

# read variables from input branch
source ./code/code.main/scripts-read-job-vars $branch "$objects" "genome genome_dir"

# run parameter script
source $params

# create path
scripts-create-path $outdir/

# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------

if ($tool == fithic) then
	# select input bedpe
	if ($bias_corrected == "TRUE")	then
        	set inputLoops = "loops_filtered_bias_cpm"
	else
        	set inputLoops = "loops_filtered_nobias_cpm"
	endif
	echo "bedpe input file is: $inputLoops"

	### annotate bedpe anchors (EE, PP, EP, PE) and generate v5cFormat file ###
	
	# filter loops
	awk 'NR>1' $branch/$object/$inputLoops.tsv | cut -f 7 > ${outdir}/qval.txt
	paste $branch/$object/$inputLoops.bedpe ${outdir}/qval.txt > ${outdir}/loops_labeled_qval.bedpe
	awk -v m=${min_anchordist} -v M=${max_anchordist} -v c=${min_activity} -v mqval=${min_qvalue} '($5-$2)>=m && ($5-$2)<=M && $7>=c && $8 <= mqval' ${outdir}/loops_labeled_qval.bedpe | cut -f 1-7 > ${outdir}/loops_labeled.bedpe
	set bedpe = ${outdir}/loops_labeled.bedpe
	set bedpe2V5C_outdir = ${outdir}/bedpe2V5C
	./code/bedpe2V5C_annot.sh ${bedpe} ${k27ac} ${tss} ${atac} ${accessible_only} ${tss_extension} ${promoter_k27ac_only} ${standarize_cpm} ${use_topLoops} ${bedpe2V5C_outdir}
	set inpfile = $outdir/bedpe2V5C/all_loops_wRev_v5cFormat.csv
else
	set inpfile = $branch/$object/virtual-5C_top200k.csv
endif

cat $inpfile | tr ',' '\t' | code/code.main/scripts-skipn 1 | awk -v D=$min_anchordist '$6>=D || $6<=-D' | sort -k8,8rg | awk -v Q="$min_qvalue" -v C="$min_activity" '$12 < Q && $8 >=C' | sort >! $outdir/loops.tsv   # apply loop filters
echo "num initital loops"
cat  $outdir/loops.tsv | wc -l
 
set main_dir = `echo ${cwd}`
cd $outdir
cat loops.tsv | grep '^ENH_' | grep 'TSS' | cut -f-2 >! ep.tsv

# calculate enhancer hubness & interactivity & promoter-associated interactivity
cat loops.tsv | grep '^ENH_' | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$1":"$3":"$4":"$5}' | sort -k 15,15 -u | cut -f1 | sort | uniq -c | tools-cols 1 0 | tr ' ' '\t' >! e.hubness.tsv
cat loops.tsv | grep '^ENH_' | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$1":"$3":"$4":"$5}' | sort -k 15,15 -u | grep 'TSS' | cut -f1 | sort | uniq -c | tools-cols 1 0 | tr ' ' '\t' >! e.p.hubness.tsv
cat loops.tsv | grep '^ENH_' | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$1":"$3":"$4":"$5}' | sort -k 15,15 -u | cut -f1,8 | sort | tools-mergeuniq -merge | tools-vectors sum -n 3 >! e.interactivity.tsv
cat loops.tsv | grep '^ENH_' | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$1":"$3":"$4":"$5}' | sort -k 15,15 -u | grep 'TSS' | cut -f1,8 | sort | tools-mergeuniq -merge | tools-vectors sum -n 3 >! e.p.interactivity.tsv

# calculate overall promoter hubness & interactivity
cat loops.tsv | grep -v '^ENH_' | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$1":"$3":"$4":"$5}' | sort -k 15,15 -u | cut -f1 | sort | uniq -c | tools-cols 1 0 | tr ' ' '\t' | sort -k1b,1 >! p.hubness.tsv
cat loops.tsv | grep -v '^ENH_' | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$1":"$3":"$4":"$5}' | sort -k 15,15 -u | cut -f1,8 | sort | tools-mergeuniq -merge | tools-vectors sum -n 3 | sort -k1b,1 >! p.interactivity.tsv

# PP and PE hubness and interactivity
cat loops.tsv | grep '^TSS_' | grep -v 'ENH_' | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$1":"$3":"$4":"$5}' | sort -k 15,15 -u | cut -f1 | sort | uniq -c | tools-cols 1 0 | tr ' ' '\t' | sort -k1b,1 >! p.p.hubness.tsv
cat loops.tsv | grep '^TSS_' | grep -v 'ENH_' | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$1":"$3":"$4":"$5}' | sort -k 15,15 -u | cut -f1,8 | sort | tools-mergeuniq -merge | tools-vectors sum -n 3 | sort -k1b,1 >! p.p.interactivity.tsv
cat loops.tsv | grep '^TSS_' | grep 'ENH_' | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$1":"$3":"$4":"$5}' | sort -k 15,15 -u | cut -f1 | sort | uniq -c | tools-cols 1 0 | tr ' ' '\t' | sort -k1b,1 >! p.e.hubness.tsv
cat loops.tsv | grep '^TSS_' | grep 'ENH_' | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$1":"$3":"$4":"$5}' | sort -k 15,15 -u | cut -f1,8 | sort | tools-mergeuniq -merge | tools-vectors sum -n 3 | sort -k1b,1 >! p.e.interactivity.tsv

# calculate enhancer-promoter contribution & exclusivity
cat loops.tsv | grep '^ENH_' | grep TSS_ | cut -f1,2,8 | sort -k2 | join -t '	'  -1 2 - p.interactivity.tsv | awk '{print $2,$1,$3/$4}' | tr ' ' '\t' >! ep.contribution.tsv
cat loops.tsv | grep '^ENH_' | grep TSS_ | cut -f1,2,8 | sort | join -t '	'  - e.p.interactivity.tsv | awk '{print $1,$2,$3/$4}' | tr ' ' '\t' >! ep.exclusivity.tsv

# calculate enhancer impact 
cat ep.contribution.tsv | cut -f1,3 | sort | tools-mergeuniq -merge | tools-vectors sum -n 3 >! e.impact.tsv

# other metrics
#cat ep.exclusivity.tsv | cut -f2- | sort | tools-mergeuniq -merge | tools-vectors sum -n 3 | sed 's/TSS_//' | sort | tools-mergeuniq -merge | tools-vectors max -n 3 | sort -k1b,1 >! p.sum_exclusivity.tsv
cat ep.exclusivity.tsv | cut -f2- | sort | tools-mergeuniq -merge | tools-vectors sum -n 3 | sort | tools-mergeuniq -merge | tools-vectors max -n 3 | sort -k1b,1 >! p.sum_exclusivity.tsv

# combine enhancer metrics and annotate by promoter
join -t '	' e.hubness.tsv e.p.hubness.tsv | join -t '	'  - e.interactivity.tsv | join -t '	' - e.p.interactivity.tsv | join -t '	'  - e.impact.tsv | join -t '	' - ep.tsv >! e.metrics.tsv

# ranked gene lists
cat e.metrics.tsv | tools-cols -t 1 6 | sed 's/TSS_//' | tools-cols -t 1 0 | sort | tools-mergeuniq -merge | tools-vectors max -n 0 | sort -k2,2rg > ! genes.ranked-by-e-hubness.tsv
cat e.metrics.tsv | tools-cols -t 2 6 | sed 's/TSS_//' | tools-cols -t 1 0 | sort | tools-mergeuniq -merge | tools-vectors max -n 0 | sort -k2,2rg > ! genes.ranked-by-e-phubness.tsv
cat e.metrics.tsv | tools-cols -t 3 6 | sed 's/TSS_//' | tools-cols -t 1 0 | sort | tools-mergeuniq -merge | tools-vectors max -n 3 | sort -k2,2rg > ! genes.ranked-by-e-interactivity.tsv
cat e.metrics.tsv | tools-cols -t 4 6 | sed 's/TSS_//' | tools-cols -t 1 0 | sort | tools-mergeuniq -merge | tools-vectors max -n 3 | sort -k2,2rg > ! genes.ranked-by-e-pinteractivity.tsv
cat e.metrics.tsv | tools-cols -t 5 6 | sed 's/TSS_//' | tools-cols -t 1 0 | sort | tools-mergeuniq -merge | tools-vectors max -n 3 | sort -k2,2rg > ! genes.ranked-by-e-impact.tsv
cat p.sum_exclusivity.tsv | sed 's/TSS_//' | sort -k2,2rg >! genes.ranked-by-e-exclusivity.tsv

awk -v OFS="\t" '{print $1":"$2,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14}' loops.tsv | sort -k1b,1 > tmp.txt; mv tmp.txt loops.tsv

# generate master EP file
awk -v OFS="\t" '{print $1":"$2,$3}' ep.contribution.tsv | sort -k1b,1 > EC.tsv
awk -v OFS="\t" '{print $1":"$2,$3}' ep.exclusivity.tsv | sort -k1b,1 > EX.tsv

#join EC.tsv EX.tsv | tr ' ' '\t' | awk -v OFS="\t" '{print $1,$2,$3,$2+$3,$2-$3}' | tr ':' '\t' | sort -k1b,1 | join - e.hubness.tsv | join - e.interactivity.tsv | join - e.impact.tsv | join - e.p.hubness.tsv | join - e.p.interactivity.tsv | sort -k2b,2 | join -1 2 - p.hubness.tsv | tr ' ' '\t' > tmp

join EC.tsv EX.tsv | tr ' ' '\t' | awk -v OFS="\t" '{print $1,$2,$3,$2+$3,$2-$3}' | tr ':' '\t' | sort -k1b,1 | join - e.hubness.tsv | join - e.interactivity.tsv | join - e.impact.tsv | join - e.p.hubness.tsv | join - e.p.interactivity.tsv | sort -k2b,2 | join -1 2 - p.hubness.tsv | join - p.p.hubness.tsv | join - p.e.hubness.tsv | join - p.p.interactivity.tsv | join - p.e.interactivity.tsv | join - p.interactivity.tsv | join - p.sum_exclusivity.tsv | tr ' ' '\t' > tmp


#cat p.interactivity.tsv | sort -k1b,1 > tmp1
#cat tmp | sort -k1b,1 > tmp2
#join tmp2 tmp1 | awk -v OFS="\t" '{print $1":"$2,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}' | sort -k1b,1 > tmp3
#join tmp2 tmp1 | awk -v OFS="\t" '{print $1":"$2,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17}' | sort -k1b,1 > tmp3

#cat loops.tsv | sort -k1b,1 > tmp4
#cat tmp3 | sort -k1b,1 > tmp5
#join tmp5 tmp4 | sort -k2b,2 > C

awk -v OFS="\t" '{print $1":"$2,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18}' tmp | sort -k1b,1 > tmp1

echo "EP.id id1 id2 EC EX EC.EX_sum EC.EX_diff e.hub e.interactivity e.impact p.EX.sum e.phub e.p.interactivity p.hub pp.hub pe.hub pp.interactivity pe.interactivity p.interactivity chr start end distance loop.activity" | tr ' ' '\t' >> EP_master.tsv

join tmp1 loops.tsv | sort -k2b,2 | awk -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$19,$11,$12,$13,$14,$15,$16,$17,$18,$22,$23,$24,$25,$27}' >> EP_master.tsv

#awk '{print "TSS_"$1"\t"$2}' p.sum_exclusivity.tsv | sort -k1b,1 > M
#join -1 2 C M | awk -v OFS='\t' '{print $2,$1,$3,$4,$5,$6,$7,$8,$9,$10,$29,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28}' > CM.txt
#awk -v OFS='\t' '{print $2,$1,$3,$4,$5,$6,$7,$8,$9,$10,$19,$11,$12,$13,$14,$15,$16,$17,$18,$22,$23,$24,$25,$27}' C > CM.txt
#cat CM.txt | tr ' ' '\t' | cut -f 1-15,18-21,23 >> EP_master.tsv
#cat CM.txt | tr ' ' '\t' | cut -f 1-19,22-25,27 >> EP_master.tsv

echo "EP.id id1 id2 EC EX EC.EX_sum EC.EX_diff e.hub e.interactivity e.impact p.EX.sum e.phub e.p.interactivity p.hub pp.hub pe.hub pp.interactivity pe.interactivity p.interactivity chr start end distance loop.activity" | tr ' ' '\t' >> EP_master_uniq.tsv
cat EP_master.tsv | sort -u -k1,1 > EP_master_uniq.tsv

# clean up
rm -f tmp* EX.tsv EC.tsv C M CM.txt

# add more metrics
$main_dir/code/calcMoreMetrics.sh $main_dir $k27ac $tss $atac

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------

# save variables
cd $main_dir
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."
