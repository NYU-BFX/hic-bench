#!/bin/tcsh

source ./inputs/params/params.tcsh

module load bwa/0.7.17

set aligner = bwa
set genome = `./code/read-sample-sheet.tcsh $sheet $object genome`
set genome_index = inputs/data/genomes/$genome/genome/bwa.index/genome.fasta
set align_params = "-A1 -B4 -E50 -L0"


