#!/bin/bash -l
#SBATCH -J hicConvert
#SBATCH --mem=40G
#SBATCH --time=12:00:00
#SBATCH -N 1
#SBATCH -c 1

cwd=$1
outdir=$2
format=$3
resolution=$4
keep_all=$5

module unload python
module unload anaconda3
module unload miniconda3
module load anaconda3/cpu/5.2.0

# Check if its a single or multiple resolution run
if [[ "$resolution" =~ "," ]]; then    multi=TRUE; fi
cd $cwd/$outdir

### COOL / MCOOL ###
if [[ $format = cool || $format = h5 || $format = homer || $format = ginteractions ]]
then
# Convert .hic to .cool/mcool
	echo 'Converting .hic to .cool file...'
	conda activate hicexplorer
	hicConvertFormat --matrices filtered.hic --outFileName filtered.cool --inputFormat hic --outputFormat cool
	conda deactivate
fi

### H5 ###
if [[ $format = h5  && $multi != TRUE || $format = ginteractions ]]
then
# Convert .cool to .h5
	echo 'Converting .cool to .h5 file...'
	conda activate hicexplorer
	hicConvertFormat --matrices filtered.cool --outFileName filtered.h5 --inputFormat cool --outputFormat h5
	conda deactivate
fi

if [[ $format = h5  && $multi = TRUE ]]
then
	echo 'hicConvertFormat doesnt support .mcool to multiple resolution .h5 format conversion'
fi

### GInteractions ###
if [[ $format = ginteractions ]]
then
# Convert .h5 to .ginteractions
        echo 'Converting .h5 to .ginteractions file...'
	hicExport -in filtered.h5 -o filtered --outputFormat GInteractions --inputFormat h5

        echo 'Generating .bedpe file...'
	cat filtered.GInteractions.tsv | tr "'" "_" | sed 's/b_/chr/g' | sed 's/_//g' | awk -v OFS="\t" '{print $1,$2+1,$3,$4,$5+1,$6,$7}' > filtered_GI.bedpe
fi

### HOMER ###
if [[ $format = homer && $multi != TRUE ]]
then
# Convert .cool to .homer
	echo 'Converting .cool to .homer file...'
	conda activate hicexplorer
	hicConvertFormat --matrices filtered.cool --outFileName filtered.homer --inputFormat cool --outputFormat homer
	conda deactivate
fi

if [[ $format = homer  && $multi = TRUE ]]
then
        echo 'hicConvertFormat doesnt support .hic to multiple resolution .homer format conversion'
fi

### Cleanup ###
rm -f filtered.reg filtered.bed 

if [[ $keep_all = FALSE && $format != cool && $format != juicer ]]
then
	rm -f filtered.hic filtered.cool
fi

if [[ $keep_all = FALSE && $format = cool ]]
then
        rm -f filtered.hic
fi

if [[ $keep_all = FALSE && $format = ginteractions ]]
then
    	rm -f filtered.h5
fi
