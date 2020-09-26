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
	# annotate bedpe anchors (EE, PP, EP, PE) and generate v5cFormat file
	set bedpe = "$branch/$object/loops_filtered_nobias_cpm.bedpe"
	./code/bedpe2V5C_annot.sh $bedpe $k27ac $tss $outdir/bedpe2V5C
	set inpfile = $outdir/bedpe2V5C/all_loops_wRev_v5cFormat.csv
else
	set inpfile = $branch/$object/virtual-5C.csv
endif

echo "num initital loops"
wc -l $inpfile

cat $inpfile | tr ',' '\t' | code/code.main/scripts-skipn 1 | awk -v D=$min_anchordist '$6>=D || $6<=-D' | sort -k8,8rg | awk -v Q="$min_qvalue" -v C="$min_activity" '$12 < Q && $8 >=C' | sort >! $outdir/loops.tsv   # apply loop filters
echo "num initital loops"
wc -l $outdir/loops.tsv
 
set main_dir = `echo ${cwd}`
cd $outdir
cat loops.tsv | grep '^ENH_' | grep 'TSS' | cut -f-2 >! ep.tsv

# calculate enhancer hubness & connectivity & promoter-associated connectivity
cat loops.tsv | grep '^ENH_' | cut -f1 | sort | uniq -c | tools-cols 1 0 | tr ' ' '\t' >! e.hubness.tsv
cat loops.tsv | grep '^ENH_' | grep 'TSS' | cut -f1 | sort | uniq -c | tools-cols 1 0 | tr ' ' '\t' >! e.phubness.tsv
cat loops.tsv | grep '^ENH_' | cut -f1,8 | sort | tools-mergeuniq -merge | tools-vectors sum -n 3 >! e.connectivity.tsv
cat loops.tsv | grep '^ENH_' | grep 'TSS' | cut -f1,8 | sort | tools-mergeuniq -merge | tools-vectors sum -n 3 >! e.pconnectivity.tsv

# calculate promoter hubness & connectivity
cat loops.tsv | grep -v '^ENH_' | cut -f1 | sort | uniq -c | tools-cols 1 0 | tr ' ' '\t' | sort -k1b,1 >! p.hubness.tsv
cat loops.tsv | grep -v '^ENH_' | cut -f1,8 | sort | tools-mergeuniq -merge | tools-vectors sum -n 3 | sort -k1b,1 >! p.connectivity.tsv

# calculate enhancer-promoter contribution & exclusivity
cat loops.tsv | grep '^ENH_' | grep TSS_ | cut -f1,2,8 | sort -k2 | join -t '	'  -1 2 - p.connectivity.tsv | awk '{print $2,$1,$3/$4}' | tr ' ' '\t' >! ep.contribution.tsv
cat loops.tsv | grep '^ENH_' | grep TSS_ | cut -f1,2,8 | sort | join -t '	'  - e.pconnectivity.tsv | awk '{print $1,$2,$3/$4}' | tr ' ' '\t' >! ep.exclusivity.tsv

# calculate enhancer impact 
cat ep.contribution.tsv | cut -f1,3 | sort | tools-mergeuniq -merge | tools-vectors sum -n 3 >! e.impact.tsv

# other metrics
cat ep.exclusivity.tsv | cut -f2- | sort | tools-mergeuniq -merge | tools-vectors sum -n 3 | sed 's/TSS_//' | sort | tools-mergeuniq -merge | tools-vectors max -n 3 | sort -k1b,1 >! p.sum_exclusivity.tsv

# combine enhancer metrics and annotate by promoter
join -t '	' e.hubness.tsv e.phubness.tsv | join -t '	'  - e.connectivity.tsv | join -t '	' - e.pconnectivity.tsv | join -t '	'  - e.impact.tsv | join -t '	' - ep.tsv >! e.metrics.tsv

# ranked gene lists
cat e.metrics.tsv | tools-cols -t 1 6 | sed 's/TSS_//' | tools-cols -t 1 0 | sort | tools-mergeuniq -merge | tools-vectors max -n 0 | sort -k2,2rg > ! genes.ranked-by-e-hubness.tsv
cat e.metrics.tsv | tools-cols -t 2 6 | sed 's/TSS_//' | tools-cols -t 1 0 | sort | tools-mergeuniq -merge | tools-vectors max -n 0 | sort -k2,2rg > ! genes.ranked-by-e-phubness.tsv
cat e.metrics.tsv | tools-cols -t 3 6 | sed 's/TSS_//' | tools-cols -t 1 0 | sort | tools-mergeuniq -merge | tools-vectors max -n 3 | sort -k2,2rg > ! genes.ranked-by-e-connectivity.tsv
cat e.metrics.tsv | tools-cols -t 4 6 | sed 's/TSS_//' | tools-cols -t 1 0 | sort | tools-mergeuniq -merge | tools-vectors max -n 3 | sort -k2,2rg > ! genes.ranked-by-e-pconnectivity.tsv
cat e.metrics.tsv | tools-cols -t 5 6 | sed 's/TSS_//' | tools-cols -t 1 0 | sort | tools-mergeuniq -merge | tools-vectors max -n 3 | sort -k2,2rg > ! genes.ranked-by-e-impact.tsv
cat p.sum_exclusivity.tsv | sed 's/TSS_//' | sort -k2,2rg >! genes.ranked-by-e-exclusivity.tsv

awk -v OFS="\t" '{print $1":"$2,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14}' loops.tsv | sort -k1b,1 > tmp.txt; mv tmp.txt loops.tsv

# generate master EP file
awk -v OFS="\t" '{print $1":"$2,$3}' ep.contribution.tsv | sort -k1b,1 > EC.tsv
awk -v OFS="\t" '{print $1":"$2,$3}' ep.exclusivity.tsv | sort -k1b,1 > EX.tsv

join EC.tsv EX.tsv | tr ' ' '\t' | awk -v OFS="\t" '{print $1,$2,$3,$2+$3,$2-$3}' | tr ':' '\t' | sort -k1b,1 | join - e.hubness.tsv | join - e.connectivity.tsv | join - e.impact.tsv | join - e.phubness.tsv | join - e.pconnectivity.tsv | sort -k2b,2 | join -1 2 - p.hubness.tsv | tr ' ' '\t' > tmp

cat p.connectivity.tsv | sort -k1b,1 > tmp1
cat tmp | sort -k1b,1 > tmp2
join tmp2 tmp1 | awk -v OFS="\t" '{print $1":"$2,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}' | sort -k1b,1 > tmp3

cat loops.tsv | sort -k1b,1 > tmp4
cat tmp3 | sort -k1b,1 > tmp5
join tmp5 tmp4 | sort -k2b,2 > C

awk '{print "TSS_"$1"\t"$2}' p.sum_exclusivity.tsv | sort -k1b,1 > M
join -1 2 C M | awk -v OFS='\t' '{print $2,$1,$3,$4,$5,$6,$7,$8,$9,$10,$29,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28}' > CM.txt

echo "EP.id id1 id2 EC EX EC.EX_sum EC.EX_diff e.hub e.interactivity e.impact p.EX.sum e.phub e.p.interactivity p.hub p.interactivity chr start end distance loop.activity" | tr ' ' '\t' >> EP_master.tsv
cat CM.txt | tr ' ' '\t' | cut -f 1-15,18-21,23 >> EP_master.tsv

echo "EP.id id1 id2 EC EX EC.EX_sum EC.EX_diff e.hub e.interactivity e.impact p.EX.sum e.phub e.p.interactivity p.hub p.interactivity chr start end distance loop.activity" | tr ' ' '\t' >> EP_master_uniq.tsv
cat EP_master.tsv | sort -u -k1,1 > EP_master_uniq.tsv
rm -f tmp* EX.tsv EC.tsv C M CM.txt

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------

# save variables
cd $main_dir
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."
