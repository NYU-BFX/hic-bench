#!/bin/tcsh

source ./inputs/params/params.tcsh

module unload gcc
module unload r
module load r/3.3.2
module load gcc/6.1.0

set chrom_excluded = 'chr[MY]'       # excluded chromosomes

set centrotelo_file = $genome_dir/centrotelo.bed

# Change these files based on your genome reference 
set gene_name = $genome_dir/gene-name_${bin_size}.tsv
set genome_file = $genome_dir/gene-info.bed

set is_normalize = 'dist_norm'
set printShowCases= FALSE

set use_sample1_ref = TRUE          # if TRUE, use SAMPLE1 TADs; if FALSE, generate common TAD boundaries
set max_boundary_dist = 3        # max distance between common boundaries (in number of bins)
set max_range = 2000000
set min_tad_size = 400000
set max_tad_size = 2e6
set max_boundary_size = 10 # max boundary length in bins.

# determine TAD branch
set tad_caller = hicratio.d_0500
set branch_short = `echo $branch | sed 's/.*results\///'`
set group_var = `echo $branch_short | cut -d'/' -f1 | cut -d'.' -f2`
set domains_branch = ../domains/results/domains.$group_var.$tad_caller/$branch_short

# Determine which data types will be integrated along with the TAD activity results
set perform_analysis = TRUE
set rnaseq = FALSE
set superenhancers = FALSE
set enhancers = FALSE
set atacseq = FALSE
