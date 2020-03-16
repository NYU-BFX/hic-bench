#!/bin/tcsh

source ./inputs/params/params.tcsh

set chrom_excluded = 'chr[MY]'                           # excluded chromosomes
set viewpoints_file = $genome_dir/protein_coding.bed     # list of coordinates to be used as viewpoints

