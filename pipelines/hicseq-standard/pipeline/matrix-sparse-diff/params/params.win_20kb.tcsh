#!/bin/tcsh

source ./inputs/params/params.tcsh

set chrom_excluded = 'chr[MY]'                           # excluded chromosomes
set viewpoints_file = $genome_dir/protein_coding.bed     # list of coordinates to be used as viewpoints

set maxdist = 2500000                                    # maximum distance from viewpoint (bp)
set window = 20000                                       # size of rolling window (bp)
set radius = `echo $window/2 | bc`                       # radius around viewpoint (bp)
set mincount = 40                                        # minimum viewpoint count for virtual 4Cs
set mindiff = 0.1                                        # minimum difference (fraction)

