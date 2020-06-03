#!/bin/bash -l
#SBATCH -J hint-tl
#SBATCH --mem-per-cpu=20G
#SBATCH --time=12:00:00
#SBATCH -N 1
#SBATCH -c 8

module unload python
module load anaconda3/cpu/5.2.0

conda activate hint
hint tl -m HICFILE --refdir /gpfs/data/skoklab/home/rodrij92/biodata/hint/hintref/hg19/ --backdir /gpfs/data/skoklab/home/rodrij92/biodata/hint/background/hg19/ --ppath pairix -f juicer -g GENOME -e ENZYME -n tl -p 8 -o ./hint/tl/
conda deactivate
