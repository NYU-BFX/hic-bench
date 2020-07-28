#!/bin/bash
#SBATCH -J APA
#SBATCH --mem=20G
#SBATCH --time=2:00:00
#SBATCH -N 1

module unload r
module load juicer/1.5

hic_file=$1
bedpes=$2
res=$3
inpdir=$4
object=$5
main_dir=$6
analysis=$7
fmin_distance=$8

bedpe=`cd $main_dir/params/; ls -l *bedpe | awk -v n="${SLURM_ARRAY_TASK_ID}" 'FNR==n {print $9}'`
bedpe_path=$main_dir/params/$bedpe
outname=`echo $bedpe | sed 's/.bedpe//g'`
outdir=$inpdir/APA/"$analysis"/"$outname"_"$object"

juicer_tools apa -r $res -n $fmin_distance -u $hic_file $bedpe_path $outdir
