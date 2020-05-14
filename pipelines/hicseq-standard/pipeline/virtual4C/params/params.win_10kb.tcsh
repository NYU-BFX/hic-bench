#!/bin/tcsh

source ./inputs/params/params.tcsh

set chrom_excluded = 'chr[MY]'                                         # excluded chromosomes
set resolution = 10000                                                 # resolution (bp)
set viewpoints_file = $genome_dir/protein-coding.bed                   # bed file: list of coordinates to be used as viewpoints
set maxdist = 2500000                                                  # maximum distance from viewpoint (bp)


