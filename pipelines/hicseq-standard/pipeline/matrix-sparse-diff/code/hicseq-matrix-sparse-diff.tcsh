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
source ./code/code.main/scripts-read-job-vars $branch "$objects" "genome genome_dir unit"

# run parameter script
source $params

# create path
scripts-create-path $outdir/

# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------

# determine number of reads in each sample
set normalization = genome-intra
scripts-send2err "normalization = $normalization"
if ($normalization == "genome-all") then
  set n_reads1 = `cat $branch/$object1/stats.tsv | grep '^read-pairs	' | cut -f2`
  set n_reads2 = `cat $branch/$object2/stats.tsv | grep '^read-pairs	' | cut -f2`
  scripts-send2err "- reads in sample 1 = $n_reads1"
  scripts-send2err "- reads in sample 2 = $n_reads2"
else if ($normalization == "genome-intra") then
  set n_reads1 = `cat $branch/$object1/stats.tsv | grep '^ds-accepted-intra	' | cut -f2`
  set n_reads2 = `cat $branch/$object2/stats.tsv | grep '^ds-accepted-intra	' | cut -f2`
  scripts-send2err "- reads in sample 1 = $n_reads1"
  scripts-send2err "- reads in sample 2 = $n_reads2"
else
  set n_reads1 = 0
  set n_reads2 = 0
endif

# determine radius around anchors
set radius = `echo $resolution/2 | bc`

# set options
set OPTIONS = "--nreads1=$n_reads1 --nreads2=$n_reads2 --unit=$unit --maxdist=$maxdist --radius=$radius --mincount=$mincount --mindiff=$mindiff"

# run analysis per chromosome
set CHR = `cat $genome_dir/genome.bed | cut -f1 | grep -wvE "$chrom_excluded"`
set jid =
foreach chr ($CHR)
  mkdir -p $outdir/$chr
  cat $viewpoints_file | gtools-regions bed | cut -f-6 | awk -v c=$chr '$1==c' >! $outdir/$chr/vp.bed         # generate chromosome-specific viewpoints file 
  cat $anchors_file | gtools-regions bed | cut -f-6 | awk -v c=$chr '$1==c' >! $outdir/$chr/anchors.bed       # generate chromosome-specific target anchors file 
  if (`cat $outdir/$chr/vp.bed | wc -l`>0) then 
    echo "Chromosome $chr..." | scripts-send2err
    set jpref = $outdir/__jdata/job.$chr
    set mem = 40G
    scripts-create-path $jpref
    set Rcmd = "Rscript ./code/sparse-matrix-diff.r $OPTIONS --vp-file=$outdir/$chr/vp.bed --target-file=$outdir/$chr/anchors.bed $outdir/$chr $chr $branch/$object1/matrix.$chr.mtx $branch/$object2/matrix.$chr.mtx $object1 $object2"
    echo $Rcmd | scripts-send2err
    set jid = ($jid `scripts-qsub-run $jpref 1 $mem $Rcmd`)
  endif
end

# wait until all jobs are completed
scripts-send2err "Waiting until all jobs are completed..."
scripts-qsub-wait "$jid"

# combine results from all chromosomes
cat $outdir/*/stats.csv | grep '^,' | sort -u | sed 's/^,/Gene,/' >! $outdir/stats.csv
cat $outdir/*/stats.csv | grep -v '^,' >> $outdir/stats.csv
set diff_files = $outdir/*/diff-regions.csv
head -1 $diff_files[1] >! $outdir/diff-regions.csv
foreach diff_file ($diff_files)
  cat $diff_file | scripts-skipn 1 >> $outdir/diff-regions.csv
end
set diff_files = $outdir/*/diff-anchors.csv
head -1 $diff_files[1] >! $outdir/diff-anchors.csv
foreach diff_file ($diff_files)
  cat $diff_file | scripts-skipn 1 >> $outdir/diff-anchors.csv
end

# remove bias and filter anchor pairs
Rscript ./code/filter-diff-anchors.r --output-dir=$outdir --min-dist=$mindist --min-val=$minval $outdir/diff-anchors.csv

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


