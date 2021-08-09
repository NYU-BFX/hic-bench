#!/bin/tcsh

source ./inputs/params/params.tcsh

set bin_size = 20000
set max_dist = 0
set ref = $genome_dir/genome.bed
set matrix_params = "--bin-size $bin_size --max-dist $max_dist -R $ref"

