#!/bin/tcsh

source ./inputs/params/params.tcsh

module unload gcc
module unload r
module load r/3.3.2
module load gcc/6.1.0

set chrom_excluded = 'chr[MY]'       # excluded chromosomes

set centrotelo_file = $genome_dir/centrotelo.bed

# Change these files based on your genome reference 
set gene_name = /gpfs/data/tsirigoslab/public/hic-bench/mm10/gene-name_${bin_size}.tsv
set genome_file = /gpfs/data/tsirigoslab/public/hic-bench/mm10/Mus_musculus.GRCm38.85_formatted.bed

set is_normalize = 'cpm'
set printShowCases= FALSE

set max_boundary_dist = 3        # max distance between common boundaries (in number of bins)
set max_range = 2000000
set min_tad_size = 400000
set max_tad_size = 2e6
set max_boundary_size = 10 # max boundary length in bins.

# Determine which data types will be integrated along with the TAD activity results
set perform_analysis = FALSE
set rnaseq = FALSE
set superenhancers = FALSE
set enhancers = FALSE
set atacseq = FALSE

# determine TAD branch
set tad_caller = hicratio.d_0500
set branch_short = `echo $branch | sed 's/.*results\///'`
set group_var = `echo $branch_short | cut -d'/' -f1 | cut -d'.' -f2`
set domains_branch = ../domains/results/domains.$group_var.$tad_caller/$branch_short


