#!/bin/tcsh

# load basic tools
module load r/3.6.1
module load python/cpu/2.7.15
module load samtools/1.3
module load bedtools/2.27.1
module load java/1.8
module load gsl/1.15
#module load gtools/3.0

# load tools required for each step of the pipeline (this can be overriden in local param scripts)
module load bowtie2/2.3.4.1

# sample sheet file
set sheet = inputs/sample-sheet.tsv

# global pipeline variables
set pipeline_max_jobs = 20

