#!/bin/tcsh

source ./inputs/params/params.tcsh

# Freely refer to 
set external_all = inputs/data_external
set external_group = inputs/data_external/group
set external_sample = inputs/data_external/sample

# Select output from domains-diff
# tad_caller: "hicratio.d_0500" and "topdom.win5" are available
set tad_caller = hicratio.d_0500
# is_normalize: "cpm" and "dist_norm" are available
set is_normalize = dist_norm
# use_sample1_ref: TRUE or FALSE. TRUE selects "ref1" and FALSE selects "common".
set use_sample1_ref = FALSE


# Note: Because of syntax in pyGenomeTracks, change the extension to ".bw" or "bigwig", not "bigWig".
# Change the entries of ${FILES[@]} variable accordingly.

set hic_matrix_files = ''
set domains_bed_files = ''
set domains_diff_bed_files = ''
set v4C_bedgraph_files = ''
set v4C_viewpoint_bed_files = ''

set atac_seq_bigwig_files = ''
set h3k27ac_chip_seq_bigwig_files = ''
set rna_seq_bigwig_files = ''

set gene_gtf_files = ''



