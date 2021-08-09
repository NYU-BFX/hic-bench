#!/bin/tcsh

source ./inputs/params/params.tcsh


set tool = fithic
set winsize = 10000                             # 10kb resolution
set chrom_excluded = 'chrM'                     # excluded chromosomes
set qval = 0.1                                 # qvalue cutoff 
set mindist = 30000                            # minimum distance cutoff
set maxdist = 10000000                         # maximum distance cutoff
set top_loops = FALSE                          # subset the top loops (cpm-ranked)
