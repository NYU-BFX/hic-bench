#!/bin/tcsh

source ./inputs/params/params.tcsh

set chrom_excluded = 'chr[MY]'                                  # excluded chromosomes
set viewpoints_file = ./params/viewpoints.bed                    # bed file: list of coordinates to be used as viewpoints
set anchors_file = ./params/anchors.bed             		# bed file: list of target anchors (e.g. enhancers)

set resolution = 5000                                           # resolution (bp)
set maxdist = 1000000                                           # maximum distance from viewpoint (bp)
set minvalue = 1.0                                              # minimum CPK2B (counts per kilobase^2 per billion reads) applied to virtual 5C results

set fdr_cut = 0.01						# binomial test cutoff value
set p_cut = 0.001						# minimum observed probability cutoff

