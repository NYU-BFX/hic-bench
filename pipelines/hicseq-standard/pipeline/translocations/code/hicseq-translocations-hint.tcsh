#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: hicseq-template.tcsh OUTPUT-DIR PARAM-SCRIPT BRANCH OBJECT(S)
##

if ($#argv != 4) then
  grep '^##' $0
  exit
endif

set outdir = $1
set params = $2
set genome = $3
set enzyme = $4

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

cp ./code/scripts-hint-tl.sh $outdir/hint-tl.sh
cd $outdir
sed 's/GENOME/'$genome'/g' hint-tl.sh > hint-tl2.sh
sed 's/ENZYME/'$enzyme'/g' hint-tl2.sh > hint-tl-custom.sh

sbatch hint-tl-custom.sh
rm -f hint-tl*sh

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------

# done
scripts-send2err "Done."
