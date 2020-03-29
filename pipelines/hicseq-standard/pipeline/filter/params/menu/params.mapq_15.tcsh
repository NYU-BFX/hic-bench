#!/bin/tcsh

source ./inputs/params/params.tcsh

set mapq = 15
set mindist = 25000
set filter_params = "--mapq $mapq --min-dist $mindist --max-offset 500"

