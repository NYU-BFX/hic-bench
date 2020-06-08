#!/bin/tcsh

module load latex/2019

### REQUIRED PARAMETERS ###
set report = hicseq-standard-report # report name
set standard = "TRUE"	# creates standard report

# set alignment & filtering parameters (you can find this info in the respective step's 'params' files)
set aligner = bowtie2	# e.g. bowtie2, bwa
set mapq = 20	
set mindist = 25000

# set matrix parameters
set binsize = 40 	# matrix resolution (kb)
set matrix = ic		# (e.g. ic, dist_norm)
set caller = hicratio 	# (e.g. hicratio, crane, topdom)
set compBinSize = 100  	# compartments resolution (kb)

# set domains-diff parameters (intraTAD activity section)
set by_group = "FALSE"	# creates intraTAD activity slides showing group pair-wise comparisons
set comparisons_group = "" # select which group-comparisons will be included" e.g. "CUTLL1_DMSO_A.CUTLL1_THZ1 CUTLL1_DMSO_H.CUTLL1_gSI"
set by_sample = "FALSE" # creates intraTAD activity slides showing sample pair-wise comparisons
set comparisons_sample = "" # select which sample-comparisons will be included

# set loops parameters
set include_loops = "FALSE"
set include_loops_diff = "FALSE"

### OPTIONAL PARAMETERS ### 

# set external data transfer parameters
set transferGdrive = "FALSE"	# transfers all the report data to a gdrive account (requires rclone-gdrive remote setup: https://rclone.org/drive/)
set gdrive_remote_path = ""	# e.g. Javi_remote:cluster2gdrive/Work/ABL-delivery/reports/
set transferExternal = "FALSE"    # only useful for ABL core member (requires special HPC permissons)
