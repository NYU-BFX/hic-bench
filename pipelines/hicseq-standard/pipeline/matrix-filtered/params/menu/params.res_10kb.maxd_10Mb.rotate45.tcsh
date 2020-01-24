#!/bin/tcsh

source ./inputs/params/params.tcsh

set bin_size = 10000
set max_dist = 10000000
set ref = $genome_dir/genome.bed
set matrix_params = "--bin-size $bin_size --max-dist $max_dist --rotate45 -R $ref"

