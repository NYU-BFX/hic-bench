#!/bin/tcsh

source ./inputs/params/params.tcsh

set chrom_excluded = 'chr[MY]'                                         # excluded chromosomes

set viewpoints_file = $genome_dir/protein-coding-tss.bed               # bed file: list of coordinates to be used as viewpoints
set anchors_file = $viewpoints_file                                    # bed file: list of target anchors (e.g. enhancers)

set resolution = 10000                                                 # resolution (bp)
set maxdist = 2500000                                                  # maximum distance from viewpoint (bp)
set mincount = 10                                                      # minimum count filter for virtual 5C


