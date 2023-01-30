#!/bin/bash
#$ -S /bin/sh
#SBATCH -J peakCall_1
#SBATCH --mem=5G
#SBATCH --time=1:00:00
#SBATCH -N 1

#### RUN EXAMPLE ########
# ./code/peakCall.sh    #
#########################

branch=`realpath results/*/*/ | fgrep -v "total"`
n_samples=`ls -l $branch | fgrep -v "total" | wc -l`
make_bam=FALSE

echo $n_samples 
echo $branch
echo $make_bam

sbatch --array=1-$n_samples ./code/scripts-peakCall.sh ${branch} ${make_bam}
