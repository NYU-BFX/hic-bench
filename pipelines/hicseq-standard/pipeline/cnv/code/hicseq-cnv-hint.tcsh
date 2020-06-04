#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: hicseq-template.tcsh OUTPUT-DIR PARAM-SCRIPT BRANCH OBJECT(S)
##

if ($#argv != 5) then
  grep '^##' $0
  exit
endif

set outdir = $1
set params = $2
set genome = $3
set enzyme = $4
set hic_file = $5

# run parameter script
source $params

# create path
scripts-create-path $outdir/

# save variables
source ./code/code.main/scripts-save-job-vars

# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------

mkdir $outdir/hint
mkdir $outdir/hint/cnv

set main_dir = `echo ${cwd}`
set hicfile_full = $main_dir/$hic_file
set refdir = $main_dir/inputs/data/genomes/$genome/hint-data/hintref/
set bicseq = $main_dir/inputs/data/tools/BICseq2-seg_v0.7.3

cp ./code/scripts-hint-cnv.sh $outdir/hint-cnv.sh

cd $outdir
sed -i 's/GENOME/'$genome'/g' hint-cnv.sh
sed -i 's/ENZYME/'$enzyme'/g' hint-cnv.sh
sed -i 's/RESOLUTION/'$resolution'/g' hint-cnv.sh
sed -i 's|REFDIR|'$refdir'|g' hint-cnv.sh                
sed -i 's|BICSEQ|'$bicseq'|g' hint-cnv.sh                
sed -i 's|HICFILE|'$hicfile_full'|g' hint-cnv.sh

echo "Calling cnv with HiNT..."
mkdir __jdata
set jid = `sbatch --output="__jdata/job.%a.out" --error="__jdata/job.%a.err" hint-cnv.sh`
set jid = `echo $jid | sed 's/.* //'`
echo $jid >! __jdata/job.id

echo "Waiting for job array [$jid] to complete..."
cd $main_dir
scripts-qsub-wait "$jid"

rm -f $outdir/hint-cnv*sh

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------

# done
scripts-send2err "Done."
