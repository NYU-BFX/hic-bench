#!/bin/tcsh

source ./inputs/params/params.tcsh

# TAD caller and intra-TAD method
set domains_method = topdom.win_5
set domains_diff_method = ${domains_method}.dist_norm.ref1

# Run template parameter script
source params/params-template.selected-regions.tcsh

