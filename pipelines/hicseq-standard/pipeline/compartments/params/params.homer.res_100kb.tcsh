#!/bin/tcsh

source ./inputs/params/params.tcsh

module load homer/4.10

set tool = homer
set resolution = 100000                               # 100kb resolution
set active_mark =           		              # active mark bed file (e.g. H3K27ac) 
set HK_genes = inputs/genomes/$genome/HK_genes.bed    # house-keeping genes bed file
set TSS_genes = inputs/genomes/$genome/gene-tss_annot.bed

