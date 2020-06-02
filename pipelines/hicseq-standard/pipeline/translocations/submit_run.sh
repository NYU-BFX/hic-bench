#!/bin/bash -l
#SBATCH -J hicbch
#SBATCH --mem=10G
#SBATCH --time=3:00:00
#SBATCH -N 1

cd /gpfs/home/rodrij92/leukemia-cell-line/translocations

./run
