#!/bin/tcsh

source ./inputs/params/params.tcsh
set tool = fithic
set k27ac = ./params/k27ac_labeled.bed
set tss = ./params/3tss_relabeled.bed

set min_activity = 0
set min_qvalue = 0.01
set min_anchordist = 0
set min_hubness = 0
