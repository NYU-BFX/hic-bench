#!/bin/tcsh

source ./inputs/params/params.tcsh
set tool = fithic
<<<<<<< HEAD
set k27ac = ./params/k27ac.bed
set tss = ./params/tss.bed
set atac = ./params/atac.bed

set accessible_only = "FALSE"
set tss_extension = 5000
set bias_corrected = "FALSE"
set promoter_k27ac_only = "FALSE"
=======
set k27ac = ./params/k27ac_labeled.bed
set tss = ./params/3tss_relabeled.bed
>>>>>>> afe0785c1090bf9222019dda27648b0b82a31f59

set min_activity = 0
set min_qvalue = 0.01
set min_anchordist = 0
set min_hubness = 0
