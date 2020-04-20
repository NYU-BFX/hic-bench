#!/bin/tcsh

source ./inputs/params/params.tcsh

set chrom_excluded = 'chr[MY]'                           # excluded chromosomes
set viewpoints_file = $genome_dir/protein_coding.bed     # bed file: list of coordinates to be used as viewpoints
set anchors_file =                                       # bed file: list of target anchors (e.g. enhancers)

set maxdist = 2500000                                    # maximum distance from viewpoint (bp)
set window = 10000                                       # size of rolling window (bp)
set radius = `echo $window/2 | bc`                       # radius around viewpoint (bp)
set mincount = 40                                        # minimum viewpoint count for virtual 4Cs
set mindiff = 0.1                                        # minimum difference (fraction)

