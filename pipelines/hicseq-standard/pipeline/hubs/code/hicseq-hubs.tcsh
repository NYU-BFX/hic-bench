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

# input loop file
set bedpe = $branch/$object/loops_filtered_nobias_cpm.bedpe

# make loops symmetric and filter by CPM and minimum distance anchor distance
( cat $bedpe ; cat $bedpe | tools-cols -t 3 4 5 0 1 2 6 ) | awk -v c=$min_cpm '$7>=c' | awk -v d=$min_dist '($2-$5>=d) || ($5-$2>=d)' | cut -f-6 | sort -u >! $outdir/filtered.bedpe
cat $outdir/filtered.bedpe | awk '{print $1","$2","$3" "$4","$5","$6}' | tools-rows -number -pref L_ | tools-key-expand | sed 's/,/ + /' | tr ',' ' ' >! $outdir/filtered.reg
echo "Found `cat $outdir/filtered.bedpe | wc -l` filtered loops." | scripts-send2err

# store all viewpoints (anchors)
cat $outdir/filtered.bedpe | tr '\t' '\n' | tools-rows -p 3 -m -t ',' | sort -u | tools-rows -number -pref VP_ | sed 's/,/ + /' | tr ',' ' ' > ! $outdir/vp.reg
echo "Found `cat $outdir/vp.reg | wc -l` viewpoints." | scripts-send2err

# generate all pairwise overlaps
cat $outdir/vp.reg | gtools-overlaps overlap -i -label $outdir/filtered.reg >! $outdir/pairwise.reg

# final hub list + info
cat $outdir/pairwise.reg | cut -f1 | sort -u | cut -d':' -f1 | uniq -cd | tools-cols 1 0 | sort >! $outdir/hubs.id
join $outdir/hubs.id $outdir/vp.reg | sed 's/ /|/' | sed 's/ /\t/' >! $outdir/hubs.reg

# classify into p-hubs
cat $outdir/hubs.reg | gtools-overlaps overlap -i -label $genome_dir/tss.bed | sed 's/:tss//' | cut -f1 | tr -s ':|' '\t' | sort -k2,2rn >! $outdir/phubs.tsv
set n_phubs = `cat $outdir/phubs.tsv | awk '$2>=5' | wc -l`
echo "Found $n_phubs p-hubs with at least 5 connections." | scripts-send2err

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------

# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."


