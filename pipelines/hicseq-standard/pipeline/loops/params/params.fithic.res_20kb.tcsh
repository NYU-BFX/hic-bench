#!/bin/tcsh

source ./inputs/params/params.tcsh

set tool = fithic
set winsize = 20000                             # 20kb resolution
set chrom_excluded = 'chrM'                     # excluded chromosomes
set qval = 0.01					# qvalue cutoff 
set mindist = 40000				# minimum distance cutoff
set maxdist = 10000000				# maximum distance cutoff

