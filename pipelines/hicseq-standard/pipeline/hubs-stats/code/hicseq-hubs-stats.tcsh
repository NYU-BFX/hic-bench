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
source ./code/code.main/scripts-read-job-vars $branch "$objects" "genome genome_dir winsize"

# run parameter script
source $params

# create path
scripts-create-path $outdir/

# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------

if ($tool == fithic) then
	# select input bedpe
	if ( $bias_corrected == "TRUE" && $cpm_normalized == "TRUE" ) then
		set inputLoops = "loops_unfiltered_bias_cpm"

	else if ( $bias_corrected == "FALSE" && $cpm_normalized == "TRUE" ) then
		set inputLoops = "loops_unfiltered_nobias_cpm"
		
	else if ( $bias_corrected == "TRUE" && $cpm_normalized == "FALSE" ) then
		set inputLoops = "loops_unfiltered_bias_raw"
	
	else if ( $bias_corrected == "FALSE" && $cpm_normalized == "FALSE" ) then
		set inputLoops = "loops_unfiltered_nobias_raw"
	else
		set inputLoops = "wrong settings!"
	endif
	echo "bedpe input file is: $inputLoops"

	### annotate bedpe anchors (EE, PP, EP, PE) and generate v5cFormat file ###
	
	## filter loops ##
	# create uncompressed version of unfiltered loops
	echo "Uncompressing unfiltered loops..." | scripts-send2err
	cat $branch/$object/$inputLoops.tsv.gz | gunzip >! ${outdir}/$inputLoops.tsv
	awk 'NR>1' ${outdir}/$inputLoops.tsv | cut -f 7 >! ${outdir}/qval.txt
	rm -f ${outdir}/$inputLoops.tsv
	paste $branch/$object/$inputLoops.bedpe ${outdir}/qval.txt >! ${outdir}/loops_labeled_qval.bedpe
	awk -v m=${min_anchordist} -v M=${max_anchordist} -v c=${min_activity} -v mqval=${min_qvalue} '($5-$2)>=m && ($5-$2)<=M && $7>=c && $8 <= mqval' ${outdir}/loops_labeled_qval.bedpe | cut -f 1-7 > ${outdir}/loops_labeled.bedpe
	set bedpe = ${outdir}/loops_labeled.bedpe
	set bedpe2V5C_outdir = ${outdir}/bedpe2V5C
	./code/bedpe2V5C_annot.sh ${bedpe} ${k27ac} ${tss} ${atac} ${accessible_only} ${tss_extension} ${promoter_k27ac_only} ${standarize_cpm} ${use_topLoops} ${winsize} ${bedpe2V5C_outdir}
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
cat loops.tsv | grep '^ENH_' | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$1":"$3":"$4":"$5}' | grep 'TSS' | sort -k 15,15 -u | cut -f1 | sort | uniq -c | tools-cols 1 0 | tr ' ' '\t' >! e.p.hubness.tsv
cat loops.tsv | grep '^ENH_' | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$1":"$3":"$4":"$5}' | sort -k 15,15 -u | cut -f1,8 | sort | tools-mergeuniq -merge | tools-vectors sum -n 3 >! e.interactivity.tsv
cat loops.tsv | grep '^ENH_' | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$1":"$3":"$4":"$5}' | grep 'TSS' | sort -k 15,15 -u | cut -f1,8 | sort | tools-mergeuniq -merge | tools-vectors sum -n 3 >! e.p.interactivity.tsv

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

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------

# clean up
rm -f qval.txt loops_labeled* loops.tsv
mv bedpe2V5C/topLoops.bedpe finalLoops.bedpe

rm -fr bedpe2V5C

# save variables
cd $main_dir
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."
