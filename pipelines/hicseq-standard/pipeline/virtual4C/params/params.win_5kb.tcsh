#!/bin/tcsh

source ./inputs/params/params.tcsh

set chrom_excluded = 'chr[MY]'                                  # excluded chromosomes
set viewpoints_file = ./params/enhancers_promoters_LOUCY.bed	# bed file: list of coordinates to be used as viewpoints
set anchors_file = ./params/enhancers_promoters_LOUCY.bed     	# bed file: list of target anchors (e.g. enhancers)

set resolution = 5000                                           # resolution (bp)
set maxdist = 1000000                                           # maximum distance from viewpoint (bp)
set minvalue = 1                                                # minimum CPK2B (counts per kilobase^2 per billion reads) applied to virtual 5C results

set nullRmvAnchors = TRUE					# remove target-anchor counts from the null distribution computation
set v4c_bdg = FALSE						# produce v4c bedgraphs for every viewpoint

