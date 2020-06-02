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

# read variables from input branch
source ./code/code.main/scripts-read-job-vars $branch "$objects" "genome genome_dir bin_size"

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

cp ./code/scripts-hint-cnv.sh $outdir/hint-cnv.sh
cd $outdir
sed 's/GENOME/'$genome'/g' hint-cnv.sh > hint-cnv2.sh
sed 's/ENZYME/'$enzyme'/g' hint-cnv2.sh > hint-cnv3.sh
sed 's/RESOLUTION/'$resolution'/g' hint-cnv3.sh > hint-cnv-custom.sh

sbatch hint-cnv-custom.sh
rm -f hint-cnv*sh

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------

# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."


