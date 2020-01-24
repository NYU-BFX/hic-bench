#!/bin/tcsh

#source /etc/profile.d/modules.csh


# unload
module unload samtools
module unload java
module unload gcc
module unload python
module unload r

# load basic tools
#module load gtools/3.0
#module load python/2.7.3
module load python/cpu/2.7.15
module load samtools/1.3
module load bedtools/2.27.1
#module load jdk
module load java/1.7
module load gcc/8.1.0
module load r/3.3.2

# load tools required for each step of the pipeline (this can be overriden in local param scripts)
module load bowtie2/2.3.4.1
#module load armatus/2014-05-19
#module load ghmm/0.9

# sample sheet file
set sheet = inputs/sample-sheet.tsv

