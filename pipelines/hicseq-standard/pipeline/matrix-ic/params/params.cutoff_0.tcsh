#!/bin/tcsh

source ./inputs/params/params.tcsh

#module unload gcc               # this is necessary in order to take care of module conflicts in our system
module unload python                # all these module commands need to be removed: mirnylib needs to be installed for latest python version
module load python/cpu/2.7.15-2

set chrom_excluded = 'chr[MY]'       # excluded chromosomes
set cutoff = 0

