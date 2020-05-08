#!/bin/tcsh

source ./inputs/params/params.tcsh

set qcut1 = 0.01		# qvalue cutoff
set qcut2 = 0.1			# qvalue cutoff for specific-common loops classification
set common_log2FC = 1		# log2FC (cpm normalized contactCounts) threshold for common loops classification
set bias_correction = TRUE	# TRUE: Uses the bias-corrected loops / FALSE: uses nobias-corrected loops
set min_distance = 40000	# set minimum distance cutoff
set max_distance = 10000000	# set maximum distance cutoff
set min_cpm = 0.01		# set minimum cpm-normalized contactCounts cutoff
