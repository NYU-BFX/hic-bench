#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: hicseq-compartments-stats.tcsh OUTDIR PARAM-SCRIPT BOUNDARY-SCORES-BRANCH [OBJECTS]
##

if (($#argv < 1) || ($#argv > 4)) then
  grep '^##' $0
  exit
endif

set outdir = $1
set params = $2
set branch = $3
set objects = ($4)

# if objects is empty, use all objects in the branch
if ("$objects" == "") set objects = `cd $branch; ls -1d *`

# read variables from input branch
source ./code/code.main/scripts-read-job-vars $branch "$objects" "genome genome_dir"

# run parameter script
source $params

# create path
scripts-create-path $outdir/


# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------

# Check for number of objects
if ($#objects < 2) then
	scripts-send2err "WARNING: More than one input objects are required. No output will be produced!"
	exit 1
endif

# Generate Compartments Plots
foreach inp_var (sample group)
  set n = `./code/read-sample-sheet2.tcsh $sheet "$objects" "$inp_var $group_var" | wc -l`
  if ($n > 0) then
    ./code/read-sample-sheet2.tcsh $sheet "$objects" "$inp_var $group_var" | awk '{print $2":"$1}' | sort -u | cut -d'-' -f$label_fields >! $outdir/labels.tsv
    break
  endif
end
echo -n "" >! $outdir/data.tsv
# echo -n "" >! $outdir/metrics.tsv
head -1 $branch/$objects[1]/pc1_metrics_summary.txt >! $outdir/metrics.tsv

foreach object ($objects)
	cat $branch/$object/compartments.scores.bedGraph | cut -f1-4 | sed '1d' | grep -vwE "$chrom_excluded" | awk '{print $1":"$2"-"$3"\t"$4}' | sed "s/\t/\t$object\t/" >> $outdir/data.tsv
	tail -n +2 $branch/$object/pc1_metrics_summary.txt >> $outdir/metrics.tsv
end

cat $outdir/data.tsv | tools-table -c -n 6 | sed 's/ *$//' | tr -s ' ' '\t' >! $outdir/pca1.matrix.tsv
Rscript ./code/scripts-compartments-stats.r -o $outdir -L $outdir/labels.tsv -m $outdir/pca1.matrix.tsv -r $outdir/metrics.tsv -c $centrotelo_file -d $delta_cut #$pca_params

rm -f $outdir/data.tsv

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------


# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."
