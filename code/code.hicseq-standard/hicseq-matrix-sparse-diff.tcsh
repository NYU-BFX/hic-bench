#!/bin/tcsh
source ./code/code.main/custom-tcshrc    # customize shell environment

##
## USAGE: hicseq-matrix-sparse-diff.tcsh OUTPUT-DIR PARAM-SCRIPT MATRIX-BRANCH OBJECT1 OBJECT2
##

if ($#argv != 5) then
  grep '^##' $0
  exit
endif

set outdir = $1
set params = $2
set branch = $3
set object1 = ($4)
set object2 = ($5)

if (($#object1>1)||($#object2>1)) then
  send2err "Error: matrix-sparse-diff operation is not implemented for multi-object grouping."
  exit 1
endif

set objects = ($object1 $object2)

# read variables from input branch
source ./code/code.main/scripts-read-job-vars $branch "$objects" "genome genome_dir"

# run parameter script
source $params

# create path
scripts-create-path $outdir/

# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------

set CHR = `cat $genome_dir/genome.bed | cut -f1 | grep -wvE "$chrom_excluded"`
set jid =
foreach chr ($CHR)
  echo "Chromosome $chr..." | scripts-send2err
  mkdir -p $outdir/$chr
  set jpref = $outdir/__jdata/job.$chr
  set mem = 40G
  scripts-create-path $jpref
  set jid = ($jid `scripts-qsub-run $jpref 1 $mem Rscript ./code/sparse-matrix-diff.r --gene-file=$viewpoints_file $outdir/$chr $chr $branch/$object1/matrix.$chr.mtx $branch/$object2/matrix.$chr.mtx $object1 $object2`)
end

# wait until all jobs are completed
scripts-send2err "Waiting until all jobs are completed..."
scripts-qsub-wait "$jid"

# combine results from all chromosomes
cat $outdir/*/stats.csv | grep '^,' | sort -u | sed 's/^,/Gene,/' >! $outdir/stats.csv
cat $outdir/*/stats.csv | grep -v '^,' >> $outdir/stats.csv

# organize virtual 4Cs into a single directory
mkdir $outdir/v4C
foreach chr ($CHR)
  mv $outdir/$chr/*-v4C.csv $outdir/v4C
  rm -rf $outdir/$chr
end

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------

# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."


