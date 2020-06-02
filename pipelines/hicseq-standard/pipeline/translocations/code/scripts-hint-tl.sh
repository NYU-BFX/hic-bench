#!/bin/bash -l
#SBATCH -J hint-tl
#SBATCH --mem=40G
#SBATCH --time=24:00:00
#SBATCH -N 1

module unload python
module load anaconda3/cpu/5.2.0

conda activate hint
hint tl -m filtered.hic --refdir /gpfs/data/skoklab/home/rodrij92/biodata/hint/hintref/hg19/ --backdir /gpfs/data/skoklab/home/rodrij92/biodata/hint/background/hg19/ --ppath pairix -f juicer -g GENOME -e ENZYME -n tl -o ./hint/tl/
conda deactivate
