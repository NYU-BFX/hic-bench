#!/bin/tcsh

source ./inputs/params/params.tcsh

# TAD caller
set domains_method = hicratio.d_0500
set domains_diff_method = ${domains_method}.dist_norm.ref1

# Run template parameter script
source params/params-template.selected-regions.tcsh

