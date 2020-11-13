#!/bin/tcsh

source ./inputs/params/params.tcsh
set tool = fithic
set k27ac = ./params/LOUCY_H3K27ac_peaks_wSignal.bed
set tss = ./params/tss_w_nmre_sorted.bed
set atac = "FALSE"

set accessible_only = "FALSE"
set tss_extension = 5000
set bias_corrected = "FALSE"
set promoter_k27ac_only = "FALSE"

set min_activity = 0
set min_qvalue = 0.01
set min_anchordist = 0
set min_hubness = 0
