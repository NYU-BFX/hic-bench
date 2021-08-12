#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: hicseq-PeakInferHiChIP-fithichip.tcsh OUTPUT-DIR HICPRO-OUTPUT GENOME BRANCH OBJECTS
##

if ($#argv != 6) then
  grep '^##' $0 
  exit   
endif
          
set outdir = $1
set params = $2
set hicpro_output = $3    
set genome = $4     # mm10
set branch = $5
set objects = ($6)


# if objects is empty, use all objects in the branch
if ("$objects" == "") set objects = `cd $branch; ls -1d *`
#
# # Run the parameter script
source $params 

#run PeakInferHiChIP.sh
set job_dir=$outdir/PeakInferHiChIP
mkdir -p $job_dir

set jid = `sbatch --output="$job_dir/job.out" --error="$job_dir/job.err" ./code/PeakInferHiChIP.sh $hicpro_output $job_dir $genome $macs` 
set jid = `echo $jid | sed 's/.* //'`
echo $jid >! $job_dir/job.id
echo "Waiting for job array [$jid] to complete..." | scripts-send2err
scripts-qsub-wait "$jid"
