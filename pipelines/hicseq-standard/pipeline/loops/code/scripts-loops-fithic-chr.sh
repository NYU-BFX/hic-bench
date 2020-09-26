#!/usr/bin/env bash
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G

set -e
shopt -s extglob

outdir="$1"
bins="$2"
winsize="$3"
chromosomes=($(echo $4))
qval="$5"
chr="${chromosomes[$SLURM_ARRAY_TASK_ID]}"
bedpe=$(realpath "$outdir"/bedpe/$chr)
bins=$(realpath "$bins"/"$chr")
main_dir=$(pwd)

echo Chromosome list = ${chromosomes[*]}
echo "Processing chromosome $chr (job id = $SLURM_ARRAY_TASK_ID)..."

mkdir -p "$outdir"/"$chr"
cd "$outdir"/"$chr"

### Conduct Fit-Hi-C ###

# Create fragment file (f1.bed)
cut -f1,2 "$bedpe" | awk -v OFS='\t' '{print $1,$2,$2}' > l1.bed
cut -f4,5 "$bedpe" | awk -v OFS='\t' '{print $1,$2,$2}'  > r1.bed
intersectBed -a "$bins" -b l1.bed -c > fl1.bed
intersectBed -a "$bins" -b r1.bed -c > fr1.bed
paste fl1.bed fr1.bed | \
    cut -f4,5,10 | \
    sed 's/_/\t/g' | \
    awk 'BEGIN {OFS="\t"} {print $1,"0",$2,$3+$4,"1"}' > f1.bed

intersectBed -a l1.bed -b "$bins" -loj > x.bed
intersectBed -a r1.bed -b "$bins" -loj > y.bed
paste x.bed y.bed | \
    cut -f7,14 | \
    sort -k1,1n -k2,2n | \
    uniq -c | \
    sed 's/^ *//' | \
    sed 's/ +\|_/\t/g' | \
    awk -v OFS='\t' '{print $2,$3,$4,$5,$1}' > i1.bed

#rm -v !("i1.bed"|"f1.bed")

# Generate gzip files for fithic
gzip i1.bed f1.bed

# Generate bias file
python3 "$main_dir"/code/HiCKRy.py -i i1.bed.gz -f f1.bed.gz -o bias.bed.gz

### Run fithic ###
module unload python/cpu/3.6.5
module load miniconda3/4.5.1	# WE HAVE TO LOAD THIS MODULE HERE TO AVOID CONFLICTS BETWEEN MINICONDA & PYTHON #
source activate fithic

# Run fithic with bias file
fithic -i i1.bed.gz -f f1.bed.gz -t bias.bed.gz -o fit_bias -r "$winsize"
zcat fit_bias/*gz | awk '(NR == 1 || $5 >= 3)' | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$8,$9}' > loops_unfiltered_bias_raw.tsv
zcat fit_bias/*gz | awk -v qval="$qval" '{if ((NR == 1) || ($5 >= 3) && ($7 <= qval)) {print}}'| awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$8,$9}' > loops_filtered_bias_raw.tsv

# Run fithic without bias file
fithic -i i1.bed.gz -f f1.bed.gz -o fit_nobias -r "$winsize"
zcat fit_nobias/*gz | awk '(NR == 1 || $5 >= 3)' | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$8,$9}' > loops_unfiltered_nobias_raw.tsv
zcat fit_nobias/*gz | awk -v qval="$qval" '{if ((NR == 1) || ($5 >= 3) && ($7 <= qval)) {print}}' | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$8,$9}' > loops_filtered_nobias_raw.tsv
