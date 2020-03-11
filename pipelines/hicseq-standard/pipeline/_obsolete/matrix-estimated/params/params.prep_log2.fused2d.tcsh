#!/bin/tcsh

source ./inputs/params/params.tcsh

set chrom_excluded = 'chr[MY]'                   # excluded chromosomes
set prep = log2
set gamma = 0.00
set min_lambda = 0.0
set max_lambda = 1.0
set n_lambda = 6
set log2_lambda = 
set algorithm = fused2D_flsa
set zone_size = `echo 2500000/$bin_size | bc`          # number of bins (max distance = 2.5Mb)
set skip_dist = 0
set split_size = 100

set hic_params = "--flsa-verbose --row-labels --skip-distance=$skip_dist --preprocess=$prep --gamma=$gamma --min-lambda=$min_lambda --max-lambda=$max_lambda --n-lambda=$n_lambda $log2_lambda --algorithm=$algorithm --zone-size=$zone_size --split-check-size=$split_size"



