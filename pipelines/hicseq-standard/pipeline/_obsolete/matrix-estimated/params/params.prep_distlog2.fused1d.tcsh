#!/bin/tcsh

source ./inputs/params/params.tcsh

set chrom_excluded = 'chr[MY]'                        # excluded chromosomes
set prep = distlog2                                   # distance-normalization for identifying specific interactions
set gamma = 0.00
set min_lambda = 0.0
set max_lambda = Inf
set n_lambda = 1
set log2_lambda = 
set algorithm = fused1D_flsa                          # use 1D approximation for lambda=Inf
set zone_size = `echo 10000000/$bin_size | bc`        # max distance = 10Mb
set skip_dist = 0
set split_size = 1e+09

set hic_params = "--flsa-verbose --row-labels --skip-distance=$skip_dist --preprocess=$prep --min-lambda=$min_lambda --max-lambda=$max_lambda --n-lambda=$n_lambda $log2_lambda --gamma=$gamma --algorithm=$algorithm --zone-size=$zone_size --split-check-size=$split_size"



