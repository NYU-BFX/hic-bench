#!/bin/tcsh
source ./code/code.main/custom-tcshrc    # customize shell environment

##
## USAGE: hicseq-virtual5C.tcsh OUTPUT-DIR PARAM-SCRIPT MATRIX-BRANCH OBJECT
##

if ($#argv != 4) then
  grep '^##' $0
  exit
endif

set outdir = $1
set params = $2
set branch = $3
set objects = ($4)
set DEBUG = false

if ($#objects>1) then
  scripts-send2err "Error: virtual4C operation is not implemented for multi-object grouping."
  exit 1
endif

set object = $objects[1]

# read variables from input branch
source ./code/code.main/scripts-read-job-vars $branch "$object" "genome genome_dir unit"

# run parameter script
source $params

# create path
scripts-create-path $outdir/

# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------

# determine scaling factor for sequencing depth
set n_reads = `cat $branch/$object/stats.tsv | grep '^ds-accepted-intra	' | cut -f2`
scripts-send2err "- number of reads = $n_reads"

# determine radius around anchors
set radius = `echo $resolution/2 | bc`

# set options
set OPTIONS = "--nreads=$n_reads --unit=$unit --maxdist=$maxdist --radius=$radius --minvalue=$minvalue"

## Check format of viewpoints/anchors files ##
cat $viewpoints_file | gtools-regions reg | gtools-regions bed | cut -f-6 >! $outdir/viewpoints.bed
cat $anchors_file | gtools-regions reg | gtools-regions bed | cut -f-6 >! $outdir/anchors.bed
foreach f (viewpoints anchors)
  if (`cat $outdir/$f.bed | gtools-regions n | awk '$2>10' | wc -l` > 0) then
    echo "Error: $f.bed file includes regions instead of single points." | scripts-send2err
    exit
  endif
  if (`cut -f4 $outdir/$f.bed | sort | uniq -d | wc -l` > 0) then
    echo "Error: $f.bed file has duplicate labels." | scripts-send2err
    exit
  endif
end

## Process V4C/V5C data per chromosome separately ##
set CHR = `cat $genome_dir/genome.bed | cut -f1 | grep -wvE "$chrom_excluded"`
set jid =
foreach chr ($CHR)
  mkdir -p $outdir/$chr
  cat $outdir/viewpoints.bed | awk -v c=$chr '$1==c' >! $outdir/$chr/vp.bed         # generate chromosome-specific viewpoints file 
  cat $outdir/anchors.bed | awk -v c=$chr '$1==c' >! $outdir/$chr/anchors.bed       # generate chromosome-specific target anchors file 
  if (`cat $outdir/$chr/vp.bed | wc -l`>0) then 
    echo "Chromosome $chr..." | scripts-send2err
    set jpref = $outdir/__jdata/job.v5c.$chr
    set mem = 20G   #`du $branch/$object/matrix.$chr.mtx | awk '{printf "%ld\n", 5+2*$1/100000}' | tools-vectors cutoff -n 0 -u -c 40`G
    scripts-create-path $jpref
    set Rcmd = "Rscript ./code/virtual4C.r $OPTIONS --vp-file=$outdir/$chr/vp.bed --target-file=$outdir/$chr/anchors.bed $outdir/$chr $chr $branch/$object/matrix.$chr.mtx $v4c_bdg"
    echo $Rcmd | scripts-send2err
    if ($DEBUG == true) then
      $Rcmd
    else
      set jid = ($jid `scripts-qsub-run $jpref 1 $mem $Rcmd`)
    endif
  endif
end

# wait until all jobs are completed
scripts-send2err "Waiting until all jobs are completed..."
scripts-qsub-wait "$jid"

if ($v4c_bdg == TRUE) then
# organize virtual 4Cs into a single directory
scripts-send2err "Organizing virtual 4Cs into a single directory..."
foreach chr ($CHR)
  if (`cat $outdir/$chr/vp.bed | wc -l`>0) then 
    mv $outdir/$chr/*.bedgraph $outdir 
    gzip $outdir/*.bedgraph
    mv $outdir/$chr/*.bedgraph $outdir
  endif
end
endif

(cd $outdir; tar cvzf bedgraphs.tgz *.bedgraph; rm -f *.bedgraph)

mkdir $outdir/bedgraph
mv $outdir/*bedgraph.gz $outdir/bedgraph/



## Compute counts per distance distribution ##
set CHR = `cat $genome_dir/genome.bed | cut -f1 | grep -wvE "$chrom_excluded"`
set jid2 =
set genome_sz = $genome_dir/genome.bed

mkdir -p $outdir/distribution

foreach chr ($CHR)
    echo "Chromosome $chr..." | scripts-send2err
    set jpref = $outdir/__jdata/job.distri.$chr
    set mem = 5G 
    scripts-create-path $jpref
    set Rcmd = "Rscript ./code/scripts-virtual5C-distCountsChr.r $outdir/distribution $chr $n_reads $unit $max_d_db $radius $random_sites $genome_sz $branch/$object/matrix.$chr.mtx"
    echo $Rcmd | scripts-send2err
    if ($DEBUG == true) then
      $Rcmd
    else
      set jid2 = ($jid2 `scripts-qsub-run $jpref 1 $mem $Rcmd`)
    endif
  endif
end


# wait until all jobs are completed
scripts-send2err "Waiting until all jobs are completed..."
scripts-qsub-wait "$jid2"

# combine distribution data files
cat $outdir/distribution/distCounts_chr* | sort -n -k1,3 -s > $outdir/distribution.tsv

# compute splines
set n_chr = `echo $CHR | tr ' ' '\n' | wc -l`
Rscript ./code/scripts-V5C-splines.r $outdir/distribution.tsv $outdir $unit $radius $min_d_db $max_d_db $resolution_db $random_sites $n_chr

### binomial test ###
scripts-send2err "Performing a binomial test for each VP/anchor connection..."
set n = `cat $outdir/chr*/virtual-5C.csv | fgrep -v "Count" | sed 's/,/\t/g' | awk -v d=$mindist '$6 > d' | wc -l`  # total number of tests: used for qvalue computation 

# Process each chromosome separately
set jid3 =
set CHR = `cat $genome_dir/genome.bed | cut -f1 | grep -wvE "$chrom_excluded"`
foreach chr ($CHR)
    if (`cat $outdir/$chr/virtual-5C.csv | wc -l`>0) then 
    echo "Chromosome $chr..." | scripts-send2err
    set jpref = $outdir/__jdata/job.binom.$chr
    set mem = 10G
    scripts-create-path $jpref
    set v5c_file = "$outdir/$chr/virtual-5C.csv"
    set Rcmd = "Rscript ./code/scripts-V5C-binomialTest.r $outdir/spline_corrected.csv $v5c_file $outdir/$chr $unit $radius $min_d_db $max_d_db $resolution_db $chr"
    echo $Rcmd | scripts-send2err
    set jid3 = ($jid3 `scripts-qsub-run $jpref 1 $mem $Rcmd`)
  endif
end

# wait until all jobs are completed
scripts-send2err "Waiting until all jobs are completed..."
scripts-qsub-wait "$jid3"

# combine results
head -n 1 $outdir/chr1/virtual-5C_btest.csv >> $outdir/virtual-5C.csv
cat $outdir/chr*/virtual-5C_btest.csv | fgrep -v "Counts" >> $outdir/virtual-5C.csv

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------

# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."
