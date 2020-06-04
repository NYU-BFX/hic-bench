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
set refdir = $main_dir/inputs/data/genomes/$genome/hint-data/hintref/
set backdir = $main_dir/inputs/data/genomes/$genome/hint-data/background/

cp ./code/scripts-hint-tl.sh $outdir/hint-tl.sh

cd $outdir
sed -i 's/GENOME/'$genome'/g' hint-tl.sh
sed -i 's/ENZYME/'$enzyme'/g' hint-tl.sh
sed -i 's|REFDIR|'$refdir'|g' hint-tl.sh              
sed -i 's|BACKDIR|'$backdir'|g' hint-tl.sh              
sed -i 's|HICFILE|'$hicfile_full'|g' hint-tl.sh

echo "Calling translocations with HiNT..."
mkdir __jdata
set jid = `sbatch --output="__jdata/job.%a.out" --error="__jdata/job.%a.err" hint-tl.sh`
set jid = `echo $jid | sed 's/.* //'`
echo $jid >! __jdata/job.id
echo "Waiting for job array [$jid] to complete..."

cd $main_dir
scripts-qsub-wait "$jid"
rm -f $outdir/hint-tl*sh

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------

# done
scripts-send2err "Done."
