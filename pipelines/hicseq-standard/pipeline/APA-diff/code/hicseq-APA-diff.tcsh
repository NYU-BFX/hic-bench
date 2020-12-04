#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: hicseq-template.tcsh OUTPUT-DIR PARAM-SCRIPT LOOP-BRANCH OBJECT1 OBJECT2
##

if ($#argv != 5) then
  grep '^##' $0
  exit
endif

set outdir = $1
set params = $2
set branch = $3
set object1 = "$4"
set object2 = "$5"
set objects = ($object1 $object2)

# read variables from input branch
source ./code/code.main/scripts-read-job-vars $branch "genome genome_dir"

# run parameter script
source $params

# create path
scripts-create-path $outdir/

# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------

# Set parameters
set main_dir = `echo ${cwd}`
set bedpes =  `cd ./params/; ls *.bedpe`
cd $main_dir
set nbed = `ls -l ./params/*.bedpe | wc -l`

### APA Analysis ###
	
# Perform analysis on the different loops subsets #

set job_dir = $outdir/__jdata_APA
mkdir -p $job_dir

set hicfile1 = inpdirs/tracks/results/tracks.by_*/*/*/"$object1"/filtered.hic
set hicfile2 = inpdirs/tracks/results/tracks.by_*/*/*/"$object2"/filtered.hic
	
echo "Computing APA scores on the loop-subsets..." | scripts-send2err
mkdir $outdir/APA
mkdir $outdir/APA/diff

set jid1 = `sbatch --array=1-$nbed --output="$job_dir/job.%a.out" --error="$job_dir/job.%a.err" ./code/APA.sh $hicfile1 "$bedpes" $APA_resolution $outdir $object1 $main_dir diff $fmin_distance`
set jid2 = `sbatch --array=1-$nbed --output="$job_dir/job.%a.out" --error="$job_dir/job.%a.err" ./code/APA.sh $hicfile2 "$bedpes" $APA_resolution $outdir $object2 $main_dir diff $fmin_distance`

set jid1 = `echo $jid1 | sed 's/.* //'`
set jid2 = `echo $jid2 | sed 's/.* //'`

echo $jid1 >! $job_dir/job1.id
echo $jid2 >! $job_dir/job2.id

scripts-send2err "Waiting for job array [$jid1] to complete..."
scripts-send2err "Waiting for job array [$jid2] to complete..."
scripts-qsub-wait "$jid1"
scripts-qsub-wait "$jid2"
rm -fr $outdir/APA/diff/*/*/*v*

echo "Performing in-house APA analysis on the loops set..." | scripts-send2err

set inbedpes = `echo $bedpes | sed 's/ /,/g'`  # comma separated list of bedpe file names
Rscript ./code/scripts-APA-diff.r $outdir/APA/diff/ $inbedpes $APA_resolution $object1 $object2 $URm

# clean up
mv $outdir/APA/diff/results_* $outdir/
rm -fr $outdir/APA/

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------

# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."
