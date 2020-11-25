#!/bin/bash -l
#SBATCH -J hicbch
#SBATCH --mem=5G
#SBATCH --time=3:00:00
#SBATCH -N 1

./run-v5c
