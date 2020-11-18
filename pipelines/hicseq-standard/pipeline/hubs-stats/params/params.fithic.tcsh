#!/bin/tcsh

source ./inputs/params/params.tcsh
set tool = fithic
set k27ac = ./params/loucy_cutll1_enhancers_wSignal.bed		# enhancer data (required)
set tss = ./params/tss_w_nmre_sorted.bed		
set atac = "FALSE"						# accessibility data (optional)

set accessible_only = "FALSE"					# use enhancers and promoters that are accessible
set tss_extension = 5000					# promoter length (upstream from TSS site)
set bias_corrected = "FALSE"					# use fithic bias corrected loops 
set promoter_k27ac_only = "FALSE"				# use promoters that have k27ac activity
set k27ac_in_TSS_anchor = "TRUE"				# use k27ac peaks that fall in anchors with TSS sites 

set min_activity = 0.2						# minimum loop contact score (cpm)
set min_qvalue = 0.01						# minimum loop qvalue
set min_anchordist = 20000				        # minimum loop distance
set max_anchordist = 10000000					# maximum loop distance
