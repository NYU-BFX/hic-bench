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

# filter loops

#cat $branch/$object/virtual-5C.csv | tr ',' '\t' | scripts-skipn 1 | sort | awk -v L=$min_loopscore '$8>=L' | awk -v D=$min_anchordist '$6>=D || $6<=-D' >! $outdir/virtual-5C.tsv
#cat $branch/$object/virtual-5C.csv | tr ',' '\t' | scripts-skipn 1 | awk -v D=$min_anchordist '$6>=D || $6<=-D' | sort -k8,8rg | head -$ntop_loops | sort >! $outdir/virtual-5C.tsv
cat $branch/$object/virtual-5C.csv | tr ',' '\t' | scripts-skipn 1 | awk -v D=$min_anchordist '$6>=D || $6<=-D' | sort -k8,8rg | awk -v Q="$min_qvalue" -v PO="$min_pobserved" '$12 < Q && $14 > PO' | sort >! $outdir/virtual-5C.tsv   # loop qvalue filtering

set main_dir = `echo ${cwd}`
cd $outdir
cat virtual-5C.tsv | grep '^ENH_' | grep '_TSS' | cut -f-2 >! ep.tsv

# calculate enhancer hubness & connectivity & promoter-associated connectivity
cat virtual-5C.tsv | grep '^ENH_' | cut -f1 | sort | uniq -c | tools-cols 1 0 | tr ' ' '\t' >! e.hubness.tsv
cat virtual-5C.tsv | grep '^ENH_' | grep '_TSS' | cut -f1 | sort | uniq -c | tools-cols 1 0 | tr ' ' '\t' >! e.phubness.tsv
cat virtual-5C.tsv | grep '^ENH_' | cut -f1,8 | sort | tools-mergeuniq -merge | tools-vectors sum -n 3 >! e.connectivity.tsv
cat virtual-5C.tsv | grep '^ENH_' | grep '_TSS' | cut -f1,8 | sort | tools-mergeuniq -merge | tools-vectors sum -n 3 >! e.pconnectivity.tsv

# calculate promoter hubness & connectivity
cat virtual-5C.tsv | grep -v '^ENH_' | cut -f1 | sort | uniq -c | tools-cols 1 0 | tr ' ' '\t' >! p.hubness.tsv
cat virtual-5C.tsv | grep -v '^ENH_' | cut -f1,8 | sort | tools-mergeuniq -merge | tools-vectors sum -n 3 >! p.connectivity.tsv

# calculate enhancer-promoter contribution & exclusivity
cat virtual-5C.tsv | grep '^ENH_' | grep _TSS_ | cut -f1,2,8 | sort -k2 | join -t '	'  -1 2 - p.connectivity.tsv | awk '{print $2,$1,$3/$4}' | tr ' ' '\t' >! ep.contribution.tsv
cat virtual-5C.tsv | grep '^ENH_' | grep _TSS_ | cut -f1,2,8 | sort | join -t '	'  - e.pconnectivity.tsv | awk '{print $1,$2,$3/$4}' | tr ' ' '\t' >! ep.exclusivity.tsv

# calculate enhancer impact 
cat ep.contribution.tsv | cut -f1,3 | sort | tools-mergeuniq -merge | tools-vectors sum -n 3 >! e.impact.tsv

# other metrics
cat ep.exclusivity.tsv | cut -f2- | sort | tools-mergeuniq -merge | tools-vectors sum -n 3 | sed 's/_TSS_[0-9]\+//' | sort | tools-mergeuniq -merge | tools-vectors max -n 3 >! p.sum_exclusivity.tsv

# combine enhancer metrics and annotate by promoter
join -t '	' e.hubness.tsv e.phubness.tsv | join -t '	'  - e.connectivity.tsv | join -t '	' - e.pconnectivity.tsv | join -t '	'  - e.impact.tsv | join -t '	' - ep.tsv >! e.metrics.tsv

# ranked gene lists
cat e.metrics.tsv | tools-cols -t 1 6 | sed 's/_TSS.*$//' | tools-cols -t 1 0 | sort | tools-mergeuniq -merge | tools-vectors max -n 0 | sort -k2,2rg > ! genes.ranked-by-e-hubness.tsv
cat e.metrics.tsv | tools-cols -t 2 6 | sed 's/_TSS.*$//' | tools-cols -t 1 0 | sort | tools-mergeuniq -merge | tools-vectors max -n 0 | sort -k2,2rg > ! genes.ranked-by-e-phubness.tsv
cat e.metrics.tsv | tools-cols -t 3 6 | sed 's/_TSS.*$//' | tools-cols -t 1 0 | sort | tools-mergeuniq -merge | tools-vectors max -n 3 | sort -k2,2rg > ! genes.ranked-by-e-connectivity.tsv
cat e.metrics.tsv | tools-cols -t 4 6 | sed 's/_TSS.*$//' | tools-cols -t 1 0 | sort | tools-mergeuniq -merge | tools-vectors max -n 3 | sort -k2,2rg > ! genes.ranked-by-e-pconnectivity.tsv
cat e.metrics.tsv | tools-cols -t 5 6 | sed 's/_TSS.*$//' | tools-cols -t 1 0 | sort | tools-mergeuniq -merge | tools-vectors max -n 3 | sort -k2,2rg > ! genes.ranked-by-e-impact.tsv
cat p.sum_exclusivity.tsv | sed 's/_TSS_[0-9]\+//' | sort -k2,2rg >! genes.ranked-by-e-exclusivity.tsv


# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------

# save variables
cd $main_dir
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."
