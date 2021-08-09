#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: hicseq-loops-fithichip.tcsh OUTPUT-DIR PARAM-SCRIPT HIC-REG-FILES GENOME BRANCH OBJECTS
##

if ($#argv != 6) then
  grep '^##' $0 
  exit   
endif
          
set outdir = $1
set params = $2 
set reg = ($3)      # allValidPairs files 
set genome = $4     # mm10
set branch = $5
set objects = ($6)


# if objects is empty, use all objects in the branch
if ("$objects" == "") set objects = `cd $branch; ls -1d *`
#
# # Run the parameter script
source $params 

#find peak file 
set peakfile = $outdir/PeakInferHiChIP/MACS2_ExtSize/out_macs2_peaks.narrowPeak

set job_dir=$outdir/FitHiChIP
mkdir -p $job_dir

#edit config template
sed -e "s|mindist|$mindist|" \
    -e "s|maxdist|$maxdist|" \
    -e "s|allValidPairs|$PWD/$reg|" \
    -e "s|peak_file|$PWD/$peakfile|" \
    -e "s|output_dir|$PWD/$job_dir|" \
    -e "s|winsize|$winsize|" \
    -e "s|q_val|$qval|" \
    -e "s|refgenome|$genome|" \
	$PWD/$fithichip_config > $outdir/config-FitHiChIP-custom


#submit job 
set jid = `sbatch --output="$job_dir/job.out" --error="$job_dir/job.err" ./code/FitHiChIP_HiCPro.sh $outdir/config-FitHiChIP-custom`
set jid = `echo $jid | sed 's/.* //'`
echo $jid >! $job_dir/job.id
echo "Waiting for job array [$jid] to complete..." | scripts-send2err
scripts-qsub-wait "$jid"
