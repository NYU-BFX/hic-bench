#!/bin/tcsh

source ./inputs/params/params.tcsh
set tool = fithic
set k27ac = ./params/k27ac_esc_wSignal.bed
set tss = ./params/3tss_relabeled.bed
set atac = ./params/atac_esc_wSignal.bed

set accessible_only = "FALSE"
set tss_extension = 5000
set bias_corrected = "FALSE"
set promoter_k27ac_only = "FALSE"
set k27ac_in_TSS_anchor = "FALSE"

set min_activity = 0
set min_qvalue = 0.01
set min_anchordist = 0
set min_hubness = 0
