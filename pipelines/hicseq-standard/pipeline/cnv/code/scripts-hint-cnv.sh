#!/bin/bash -l
#SBATCH -J hint-cnv
#SBATCH --mem=10G
#SBATCH --time=12:00:00
#SBATCH -N 1

module unload python
module load anaconda3/cpu/5.2.0

conda activate hint
hint cnv -m HICFILE -f juicer --refdir /gpfs/data/skoklab/home/rodrij92/biodata/hint/hintref/GENOME/ -r RESOLUTION -g GENOME -n cnv -o ./hint/cnv/ --bicseq /gpfs/home/rodrij92/home/install/BICseq2-seg_v0.7.3 -e ENZYME
conda deactivate
