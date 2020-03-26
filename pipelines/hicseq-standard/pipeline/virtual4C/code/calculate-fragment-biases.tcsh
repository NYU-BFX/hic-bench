#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: calculate-fragment-biases.tcsh OUTPUT-DIR FRAGMENTS-BED GENOME-DIR WINDOW
##

if ($#argv != 4) then
  grep '^##' $0
  exit
endif

set outdir = $1          # Arima
set frag = $2            # genome/Arima.fragments.bed
set genome_dir = $3      # genome
set window = $4          # 20000

rm -rf $outdir
mkdir -p $outdir

# Analyze one chromosome at a time
set CHR = `cat $genome_dir/genome.bed | cut -f1`
foreach chr ($CHR)
  set outmat = $outdir/matrix.$chr.mtx
  echo "Storing fragment coordinates in $chr at 100bp resolution..." | scripts-send2err

  # Obtain matrix size (in 100bp resolution; 100bp is hardcoded)
  set n = `cat $genome_dir/genome.bed | awk -v c=$chr '{if ($1==c) print $3+100}' | sed 's/..$//'`
  
  # Obtain matrix data
  cat $frag | awk -v c=$chr '$1==c' | cut -f2,3 | tr '\t' '\n' | awk '$1>100' | sed 's/..$//' | awk -v n=$n '$1<=n' | sort | uniq -c | sed 's/^ *//' >! $outmat.out
  set N = `cat $outmat.out | wc -l`

  # Generate sparse matrices
  echo '%%MatrixMarket matrix coordinate integer general' >! $outmat
  echo 1 $n $N >> $outmat
  cat $outmat.out | awk '{print 1,$2,$1}' >> $outmat
  rm -f $outmat.out

  # Generate bedgraph biases file
  echo "Generating bedgraph biases file using a sliding window of ${window}bp..." | scripts-send2err
  Rscript ./code/calculate-fragment-biases.r -w $window $outdir/$chr.bedgraph $outmat
  rm -f $outmat
  
end

echo "Combining chromosomes into a single bedgraph (printing 1 every 10 entries only)..." | scripts-send2err
(echo "track type=bedGraph name=$outdir"; cat $outdir/*.bedgraph | grep -v '^track' | awk 'NR % 10 == 0') | gzip >! $outdir/biases.bedgraph.gz
rm -f $outdir/*.bedgraph 

