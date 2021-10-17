#!/bin/tcsh

source ./inputs/params/params.tcsh


set aligner = hicpro
set genome = `./code/read-sample-sheet.tcsh $sheet $object genome`
set enzyme = `./code/read-sample-sheet.tcsh $sheet $object enzyme`

set genome_index = inputs/genomes/$genome/bowtie2.index #hicpro can't find this if start with inputs/
set enzyme_path = inputs/genomes/$genome/${enzyme}.fragments.bed
set align_params_global = "--very-sensitive -L 30 --score-min L,-0.6,-0.2 --end-to-end --reorder"
set align_params_local = "--very-sensitive -L 20 --score-min L,-0.6,-0.2 --end-to-end --reorder"
set hicpro_config = inputs/config-hicpro.txt
