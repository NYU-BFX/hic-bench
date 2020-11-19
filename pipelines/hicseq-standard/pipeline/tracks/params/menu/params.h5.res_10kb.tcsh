#!/bin/tcsh

source ./inputs/params/params.tcsh

set format = h5				# juicer, cool, h5
set resolution = 10000			# comma separated list of resolutions, or default (2.5M, 1M, 500K, 250K, 100K, 50K, 25K, 10K, and 5K)
set keep_all = FALSE	                # keep intermediate matrices (juicer (hic) -> cool -> h5)
