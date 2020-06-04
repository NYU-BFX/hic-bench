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

cp ./code/scripts-hint-cnv.sh $outdir/hint-cnv.sh

cd $outdir
sed 's/GENOME/'$genome'/g' hint-cnv.sh > hint-cnv2.sh
sed 's/ENZYME/'$enzyme'/g' hint-cnv2.sh > hint-cnv3.sh
sed 's/RESOLUTION/'$resolution'/g' hint-cnv3.sh > hint-cnv4.sh
sed 's|HICFILE|'$hicfile_full'|g' hint-cnv4.sh > hint-cnv-custom.sh

echo "Calling cnv with HiNT..."
mkdir __jdata
set jid = `sbatch --output="__jdata/job.%a.out" --error="__jdata/job.%a.err" hint-cnv-custom.sh`
set jid = `echo $jid | sed 's/.* //'`
echo $jid >! __jdata/job.id

echo "Waiting for job array [$jid] to complete..."
cd $main_dir
scripts-qsub-wait "$jid"

#rm -f $outdir/hint-cnv*sh

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------

# done
scripts-send2err "Done."
