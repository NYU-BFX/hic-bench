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

# read variables from input branch
source ./code/code.main/scripts-read-job-vars $branch "$objects" "genome genome_dir bin_size"

# run parameter script
source $params

# create path
scripts-create-path $outdir/

# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------

# Error checking
if ($#objects > 1) then
  echo "ERROR: this operation allows only 1 input object." | scripts-send2err
  exit 1
endif

# Identify branches and objects that will be used as inputs
set object = $objects[1]
set sample1 = `echo $object | cut -d'.' -f1`
set sample2 = `echo $object | cut -d'.' -f2`
set branch_short = `echo $branch | sed 's/.*results\///' | sed 's/^domains-diff.[^/]\+\///'`
set matrix_step = `echo $branch_short | cut -d'/' -f1 | cut -d'.' -f1`
set matrix_branch = ../$matrix_step/results/$branch_short

# regions to plot
cat $branch/$object/final_results.tsv | grep -v logFC | awk -v fdr=$ddiff_fdr '$11<fdr' | awk -v l2fc=$ddiff_l2fc '$9>l2fc || $9<-l2fc' | cut -f1,3,5 | tr '\t' ' '  | sed 's/ /\t/' | tools-vectors format -n 0 | sed 's/ /\t/' | gtools-regions bounds -g $genome_dir/genome.bed >! $outdir/selected_regions.bed
set n = `cat $outdir/selected_regions.bed | wc -l`
echo "Selected $n regions for plotting." | scripts-send2err

set k = 1
while ($k <= $n) 
  # Pick k-th region
  set region = `cat $outdir/selected_regions.bed | head -$k | tail -1`

  # Determine all inputs
  set chr = `echo $region | cut -d' ' -f1`
  set matrix_top = "$matrix_branch/$sample1/matrix.$chr.tsv"
  set matrix_bottom = "$matrix_branch/$sample2/matrix.$chr.tsv"
  set chr_start = `echo $region | cut -d' ' -f2`
  set chr_end = `echo $region | cut -d' ' -f3`
  set sample_top = $sample1
  set sample_bottom = $sample2
  set outname = ${chr}_${chr_start}_${chr_end}

  # Print diagnostics
  echo $outname
  echo -n "${sample_top}: "
  ls -1 $matrix_top
  echo -n "${sample_bottom}: "
  ls -1 $matrix_bottom

  # Plot
  echo "Plotting..."
  Rscript ./code/differential_rect_binding_for_locus.R --normalize --binsize $bin_size $matrix_top $matrix_bottom $chr_start $chr_end $sample_top $sample_bottom $outdir/$outname

  @ k ++
end

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------

# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."


