#!/bin/tcsh

source ./inputs/params/params.tcsh

#set bedpe = "/gpfs/home/rodrij92/leukemia-hic/pipeline/loops/results/loops.by_group.fithic.res_20kb/filter.by_sample.mapq_15/align.by_sample.bowtie2/CUTLL1_DMSO_A/loops_filtered_bias_cpm.bedpe"
set APA_resolution = 10000		     # set APA resolution
set URm = 2.5			             # Upper Right Multiplier: to set the maximum value in the heatmap color scale [ UR = URm * mean(upper right quadrant counts) ]    
