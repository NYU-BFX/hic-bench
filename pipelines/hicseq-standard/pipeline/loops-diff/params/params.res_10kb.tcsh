#!/bin/tcsh

source ./inputs/params/params.tcsh

set top_loops = FALSE		# use top loops (Options: integer or FALSE)
set cpmScaled = FALSE		# scale cpm values so the mean cpm value of the sample is 1 (Options: TRUE or FALSE)
set qcut1 = 0.01		# qvalue cutoff
set qcut2 = 0.1			# qvalue cutoff for specific-common loops classification
set common_log2FC = 1		# log2FC (cpm normalized contactCounts) threshold for common loops classification
set bias_correction = FALSE	# TRUE: Uses the bias-corrected loops / FALSE: uses nobias-corrected loops
set min_distance = 80000	# set minimum distance cutoff
set max_distance = 10000000	# set maximum distance cutoff
set min_cpm = 0			# set minimum cpm-normalized contactCounts cutoff
set APA_diff = FALSE		# run in-house APA analysis on the specific-common loop subsets
set APA_quantiles = FALSE	# run in-house APA analysis by distance-quantiles (ref1)
set APA_qfile = common		# input loops file used for the quantile analysis (options: ref1 ; common)
set APA_resolution = 10000	# set APA resolution
set URm = 2.5			# Upper Right Multiplier: to set the maximum value in the heatmap color scale [ UR = URm * mean(upper right quadrant counts) ]
