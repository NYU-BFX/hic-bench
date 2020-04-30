#!/usr/bin/bash

outdir=$1
branch=$2
mat=$3
object1=$4
object2=$5

inpdir1=$branch/$object1
inpdir2=$branch/$object2

module unload python/cpu
module unload anaconda3/cpu
module unload miniconda3

module load anaconda3/cpu/5.2.0
conda activate tadsplimer

echo "outdir is: $outdir"
echo "branch is: $branch"
echo "object1 is: $object1"
echo "object2 is: $object2"
echo "mat is: $mat"

>&2 echo "Before loading TADsplimer's R, it is observed in the \$PATH of curren
t working shell that:"
>&2 echo "the bash variable R_LIBS_SITE is $R_LIBS_SITE"
if [ -s $R_LIBS_SITE ]
then
  >&2 echo "TADsplimer in conda environment has a peculiar R version."
  >&2 echo "It's best to use conda environment R's internal packages."
  R_LIBS_SITE=""
fi


Rscript ./code/TADsplimer-prep.R $mat $inpdir1 $object1 $outdir "first"
Rscript ./code/TADsplimer-prep.R $mat $inpdir2 $object2 $outdir "second"
python ./code/TADsplimer.py split_TADs -c $outdir/$object1/first_$mat,$outdir/$object2/second_$mat --contact_maps_aliases $object1,$object2 -o $outdir/$mat

module unload anaconda3/cpu

# python src/TADsplimer.py split_TADs -c example/simulation_merge.txt,example/simulation_split.txt --contact_maps_aliases first,second -o output

