#!/bin/tcsh

source ./inputs/params/params.tcsh

set tool = virtual5C

# set input data
set k27ac = ./params/k27ac_peaks.bed				# h3k27ac data (required)
set tss = ./params/tss.bed					# tss data (required)
set atac = "FALSE"						# accessibility data (optional)

# set hub parameters
set accessible_only = "FALSE"					# use enhancers and promoters that are accessible
set tss_extension = 5000					# promoter length (upstream from TSS site)
set promoter_k27ac_only = "FALSE"				# use promoters that have k27ac activity
set k27ac_in_TSS_anchor = "TRUE"				# use k27ac peaks that fall in anchors with TSS sites 

# set loop parameters
set min_activity = 0						# minimum loop contact score (cpm)
set min_qvalue = 0.1						# minimum loop qvalue
set bias_corrected = "FALSE"                                    # use fithic bias corrected loops
set cpm_normalized = "TRUE"
set min_anchordist = 30000				        # minimum loop distance
set max_anchordist = 25000000					# maximum loop distance
set standarize_cpm = "FALSE"					# standarize cpm values so the loop mean cpm value is 1
set use_topLoops = "FALSE"					# only use the selected top loops (cpm ranked)
