#!/bin/tcsh

source ./inputs/params/params.tcsh

set chrom_excluded = 'chr[MY]'                           # excluded chromosomes
set viewpoints_file = $genome_dir/protein_coding.bed     # list of coordinates to be used as viewpoints
set radius = 5000                                        # radius around viewpoint (bp)
set win = 5000                                           # size of rolling window for Hi-C count aggregation (bp)
set maxdist = 2500000                                    # maximum distance from viewpoint (bp)


