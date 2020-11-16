#!/bin/tcsh

source ./inputs/params/params.tcsh

set chrom_excluded = 'chr[MY]'                                  # excluded chromosomes
set viewpoints_file = ./params/tss_k27ac.bed     # bed file: list of coordinates to be used as viewpoints
set anchors_file = ./params/tss_k27ac.bed 		    	# bed file: list of target anchors (e.g. enhancers)

set resolution = 5000                                           # resolution (bp)
set maxdist = 1000000                                           # maximum distance from viewpoint (bp)
set minvalue = 1                                                # minimum CPK2B (counts per kilobase^2 per billion reads) applied to virtual 5C results

set resolution_db = 5000
set min_d_db = 20000
set max_d_db = 1000000
set random_sites = 50

set v4c_bdg = FALSE						# produce v4c bedgraphs for every viewpoint
set mindist = 20000						# minimum distance filter for the binomial test
set qval_cut = 0.01						# qvalue filter for the binomial test
