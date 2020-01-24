#!/usr/bin/tcsh
 
module unload perl
module load perl/5.28.0
 
set tool = crane
set cranepath = ./code/crane
set inssqr = 500000   # insulation square
set idspan = 200000   # insulation delta span
set insmode = mean    # insulation mode
set noise_thr = 0.1   # noise threshold
set bmoerr = 3        # boundary margin of error
set chrom_excluded = 'chr[MY]' # excluded chromosomes
 
set hicmatrix_params = "--min-lambda=0.0 --max-lambda=5.0 --n-lambda=6 --gamma=0"       # only applies if RData is estimated using --max-lambda=Inf
