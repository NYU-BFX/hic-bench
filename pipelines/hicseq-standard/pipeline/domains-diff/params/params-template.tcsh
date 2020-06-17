#!/bin/tcsh

source ./inputs/params/params.tcsh

set chrom_excluded = 'chr[MY]'       # excluded chromosomes

set centrotelo_file = $genome_dir/centrotelo.bed
set gene_tss = $genome_dir/gene-tss_annot.bed

# basic params
set printShowCases= FALSE
set max_boundary_dist = 3        # max distance between common boundaries (in number of bins)
set max_range = 2000000
set min_tad_size = 400000
set max_tad_size = 2e6
set max_boundary_size = 10 # max boundary length in bins.
set activity_lfc_cutoff = 0.25
set activity_fdr_cutoff	= 0.1
set activity_mdiff_cutoff = 0.1

# determine TAD branch
set branch_short = `echo $branch | sed 's/.*results\///'`
set group_var = `echo $branch_short | cut -d'/' -f1 | cut -d'.' -f2`
set domains_branch = ../domains/results/domains.$group_var.$tad_caller/$branch_short

# Determine which data types will be integrated along with the TAD activity results
set perform_analysis = TRUE
set rnaseq = FALSE
set superenhancers = FALSE
set enhancers = FALSE
set atacseq = FALSE

