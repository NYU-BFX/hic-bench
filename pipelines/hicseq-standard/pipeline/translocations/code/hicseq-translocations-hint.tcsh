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
mkdir $outdir/hint/tl
set main_dir = `echo ${cwd}`
set hicfile_full = $main_dir/$hic_file

cp ./code/scripts-hint-tl.sh $outdir/hint-tl.sh

cd $outdir
sed 's/GENOME/'$genome'/g' hint-tl.sh > hint-tl2.sh
sed 's/ENZYME/'$enzyme'/g' hint-tl2.sh > hint-tl3.sh
sed 's|HICFILE|'$hicfile_full'|g' hint-tl3.sh > hint-tl-custom.sh

echo "Calling translocations with HiNT..."
mkdir __jdata
set jid = `sbatch --output="__jdata/job.%a.out" --error="__jdata/job.%a.err" hint-tl-custom.sh`
set jid = `echo $jid | sed 's/.* //'`
echo $jid >! __jdata/job.id
echo "Waiting for job array [$jid] to complete..."

rm -f hint-tl*sh

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------

# done
cd $main_dir
scripts-send2err "Done."
