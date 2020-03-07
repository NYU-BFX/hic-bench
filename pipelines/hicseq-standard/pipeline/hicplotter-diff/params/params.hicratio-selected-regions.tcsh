#!/bin/tcsh

source ./inputs/params/params.tcsh

# TAD caller
set domains_method = hicratio.d_0500

# Domains-diff method and cutoffs
set domains_diff_method = ${domains_method}.dist_norm.ref1
set ddiff_fdr = 0.001
set ddiff_l2fc = 0.3

# Run template parameter script
source params/params-template.selected-regions.tcsh

