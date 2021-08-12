#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: hicseq-filter-stats.tcsh OUTPUT-DIR PARAM-SCRIPT FILTER-BRANCH [OBJECTS]
##

# process command-line inputs
if (($#argv < 3) || ($#argv > 4)) then
  grep '^##' $0 | scripts-send2err
  exit 1
endif

set outdir = $1
set params = $2
set branch = $3
set objects = ($4)

# if samples is empty, use all objects in the branch
if ("$objects" == "") set objects = `cd $branch; ls -1d *`
set sample_paths = `echo "$branch\t$objects" | tools-key-expand | tr '\t' '/'`

# run parameter script
source $params

# create path
scripts-create-path $outdir/

# Plot barplots
Rscript ./code/hicseq-filter-stats.r $outdir "$sample_paths"

#symbolic link to hicpro output 
if (`echo $branch | cut -f 5 -d"/"` == 'align.by_sample.hicpro') then 
   
  foreach f (`cd $branch; ls -1d *`) 
     ln -s ../../../../../inpdirs/filter/inpdirs/align/results/align.by_sample.hicpro/$f/hic_results/pic/$f/ $outdir
  end
 
endif 

# save variables
set >! $outdir/job.vars.tsv

scripts-send2err "Done."

