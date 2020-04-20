#!/bin/tcsh

source ./inputs/params/params.tcsh

set chrom_excluded = 'chr[MY]'                           # excluded chromosomes
set viewpoints_file = $genome_dir/protein_coding.bed     # list of coordinates to be used as viewpoints
set win = 5000                                           # size of rolling window for Hi-C count aggregation (bp)
set radius = `echo $win/2 | bc`                          # radius around viewpoint (bp)
set maxdist = 2500000                                    # maximum distance from viewpoint (bp)


