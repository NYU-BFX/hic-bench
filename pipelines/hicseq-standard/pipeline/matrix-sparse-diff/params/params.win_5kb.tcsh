#!/bin/tcsh

source ./inputs/params/params.tcsh

set chrom_excluded = 'chr[MY]'                           # excluded chromosomes

set viewpoints_file = $genome_dir/protein-coding.bed     # bed file: list of coordinates to be used as viewpoints
set anchors_file = $viewpoints_file                      # bed file: list of target anchors (e.g. enhancers)

set resolution = 5000                                    # resolution (bp)
set maxdist = 2500000                                    # maximum distance from viewpoint (bp)
set mindist = 10000                                      # minimum distance used to filter anchor pairs (bp)
set minval = 4.0                                         # minimum value used to filter anchor pairs (CPK2B)
set mincount = 40                                        # minimum viewpoint count for virtual 4Cs
set mindiff = 0.1                                        # minimum difference (fraction)

