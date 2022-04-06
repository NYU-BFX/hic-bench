#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: hicseq-tracks-h5.tcsh OUTPUT-DIR HIC-REG-FILES GENOME PARAMS
##
## Example: ./hicseq-tracks-h5.tcsh mESC_J1-HindIII "mESC_J1-HindIII-rep1*/filtered_reads.reg+" mm10 $params
##

if ($#argv != 4) then
  grep '^##' $0
  exit
endif

set outdir = $1
set reg = ($2)      # *.reg.gz files
set genome = $3     # e.g. mm10
set params = $4

# run parameter script
source $params

# Checks if the resolution selected is the  default. If so, it uses the default resolution values in Juicer.
if ($resolution == default) then
	scripts-send2err "Default resolution set selected: 5000,10000,25000,50000,100000,250000,500000,1000000,2500000."
	set resolution = 5000,10000,25000,50000,100000,250000,500000,1000000,2500000
else
	  scripts-send2err "Resolution/s: $resolution."
endif

# Create uncompressed version of mapped read pairs
echo 'Uncompressing filtered pairs...' | scripts-send2err
cat $reg | gunzip >! $outdir/filtered.reg

#shuf -n 5000000 $outdir/filtered.reg >! $outdir/filtered2.reg
#rm -f $outdir/filtered.reg
#mv $outdir/filtered2.reg $outdir/filtered.reg

echo 'Creating filtered.bed file...' | scripts-send2err
# Convert to bed format
awk ' BEGIN { OFS="\t"; strand["-"]="1"; strand["+"]="0" } {                    \
    if ($2 < $6)                                                                \
        print $1, strand[$3], $2, $4, 0, strand[$7], $6, $8, 1, 0, 1;           \
    else                                                                        \
        print $1, strand[$7], $6, $8, 1, strand[$3], $2, $4, 0, 0, 1;           \
}' $outdir/filtered.reg | sort -k3,3d -k7,7d >! $outdir/filtered.bed

# Generate .hic file
echo 'Generating .hic file...' | scripts-send2err
module unload r
module load juicer/1.5
java -Xmx8g -jar /gpfs/share/apps/juicer/1.5/scripts/juicer_tools.jar pre $outdir/filtered.bed $outdir/filtered.hic "$genome" -r $resolution -n
module unload juicer/1.5

# Normalize .hic file (KR)
echo 'Normalizing .hic file...' | scripts-send2err
module unload r
module load juicer/1.5
java -Xmx8g -jar /gpfs/share/apps/juicer/1.5/scripts/juicer_tools.jar addNorm -k KR $outdir/filtered.hic
module unload juicer/1.5

# Generate .cool and .h5 files (had to separate this step in a second bash script because the conda module doesn't work well on tcsh)
if ($format == cool || $format == h5 || $format == homer || $format == ginteractions) then
	echo 'Running hicConvertFormat...' | scripts-send2err
	set job_dir = $outdir/__jdata
	mkdir -p $job_dir
	set jid = `sbatch --output="$job_dir/job_hicConvertFormat.out" --error="$job_dir/job_hicConvertFormat.err" ./code/scripts-tracks-hicConvertFormat.sh $cwd $outdir $format $resolution $keep_all`
	set jid = `echo $jid | sed 's/.* //'`
endif

if ($keep_all == FALSE && $format == juicer) then
    	rm -f $outdir/filtered.reg $outdir/filtered.bed
endif
