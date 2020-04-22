#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: hicseq-template.tcsh OUTPUT-DIR PARAM-SCRIPT LOOP-BRANCH OBJECT1 OBJECT2
##

if ($#argv != 5) then
  grep '^##' $0
  exit
endif

set outdir = $1
set params = $2
set branch = $3
set object1 = $4
set object2 = $5

set objects = ($object1 $object2)

# read variables from input branch
source ./code/code.main/scripts-read-job-vars $branch "$objects" "genome genome_dir winsize"

# run parameter script
source $params

# create path
scripts-create-path $outdir/

# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------

set loops1 = $branch/$object1/loops_filtered_bias_cpm.tsv
set loops2 = $branch/$object2/loops_filtered_bias_cpm.tsv
set tss = $genome_dir/tss.bed
Rscript ./code/scripts-loops-diff.r $loops1 $loops2 $object1 $object2 $outdir $winsize $tss $common_log2FC

mkdir -p $outdir/bedpe_files $outdir/igv_files

# create bedpe files
awk -v var="$winsize" '{ if ((NR>1)) print $1"\t"($2-var/2)"\t"($2+var/2)"\t"$3"\t"($4-var/2)"\t"($4+var/2)"\t"$5}' $outdir/"$object1"_specific_loops.tsv >! $outdir/bedpe_files/"$object1"_specific_loops.bedpe
awk -v var="$winsize" '{ if ((NR>1)) print $1"\t"($2-var/2)"\t"($2+var/2)"\t"$3"\t"($4-var/2)"\t"($4+var/2)"\t"$5}' $outdir/"$object2"_specific_loops.tsv >! $outdir/bedpe_files/"$object2"_specific_loops.bedpe
awk -v var="$winsize" '{ if ((NR>1)) print $1"\t"($2-var/2)"\t"($2+var/2)"\t"$3"\t"($4-var/2)"\t"($4+var/2)"\t"$5}' $outdir/common.loops_decreased.tsv >! $outdir/bedpe_files/common.loops_decreased.bedpe
awk -v var="$winsize" '{ if ((NR>1)) print $1"\t"($2-var/2)"\t"($2+var/2)"\t"$3"\t"($4-var/2)"\t"($4+var/2)"\t"$5}' $outdir/common.loops_increased.tsv >! $outdir/bedpe_files/common.loops_increased.bedpe
awk -v var="$winsize" '{ if ((NR>1)) print $1"\t"($2-var/2)"\t"($2+var/2)"\t"$3"\t"($4-var/2)"\t"($4+var/2)"\t"$5}' $outdir/common.loops.tsv >! $outdir/bedpe_files/common.loops.bedpe

# create igv files
awk '{if(NR>1) print $1"\t"$2"\t"$4"\t\.\t1\.0"}' $outdir/common.loops_decreased.tsv | sed -e '1itrack graphType=junctions' | sort -k2 -n >! $outdir/igv_files/common.loops_decreased.igv.bed
awk '{if(NR>1) print $1"\t"$2"\t"$4"\t\.\t1\.0"}' $outdir/common.loops_increased.tsv | sed -e '1itrack graphType=junctions' | sort -k2 -n >! $outdir/igv_files/common.loops_increased.igv.bed
awk '{if(NR>1) print $1"\t"$2"\t"$4"\t\.\t1\.0"}' $outdir/common.loops.tsv | sed -e '1itrack graphType=junctions' | sort -k2 -n > ! $outdir/igv_files/common.loops.igv.bed
awk '{if(NR>1) print $1"\t"$2"\t"$4"\t\.\t1\.0"}' $outdir/"$object1"_specific_loops.tsv | sed -e '1itrack graphType=junctions' | sort -k2 -n >! $outdir/igv_files/"$object1"_specific_loops.igv.bed
awk '{if(NR>1) print $1"\t"$2"\t"$4"\t\.\t1\.0"}' $outdir/"$object2"_specific_loops.tsv | sed -e '1itrack graphType=junctions' | sort -k2 -n >! $outdir/igv_files/"$object2"_specific_loops.igv.bed

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------

# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."
