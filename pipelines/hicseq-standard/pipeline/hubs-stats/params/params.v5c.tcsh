#!/bin/tcsh

source ./inputs/params/params.tcsh
set tool = virtual5C
set k27ac = ./params/k27ac.bed					# enhancer data (required)
set tss = ./params/tss.bed		
set atac = "FALSE"						# accessibility data (optional)

set accessible_only = "FALSE"					# use enhancers and promoters that are accessible
set tss_extension = 5000					# promoter length (upstream from TSS site)
set bias_corrected = "FALSE"					# use fithic bias corrected loops 
set promoter_k27ac_only = "FALSE"				# use promoters that have k27ac activity
set k27ac_in_TSS_anchor = "TRUE"				# use k27ac peaks that fall in anchors with TSS sites 

set min_activity = 0						# minimum loop contact score (cpm)
set min_qvalue = 0.1						# minimum loop qvalue
set min_anchordist = 20000				        # minimum loop distance
set max_anchordist = 10000000					# maximum loop distance
