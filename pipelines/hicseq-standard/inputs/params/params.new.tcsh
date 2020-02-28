#!/bin/tcsh

# load basic tools
module load r/3.6.1
module load python/cpu/3.6.5
module load samtools/1.9
module load bedtools/2.27.1
module load java/1.8
module load gsl/2.5
#module load gtools/3.0

# load tools required for each step of the pipeline (this can be overriden in local param scripts)
module load bowtie2/2.3.5.1

# sample sheet file
set sheet = inputs/sample-sheet.tsv

