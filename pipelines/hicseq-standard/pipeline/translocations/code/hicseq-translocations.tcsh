#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: hicseq-translocations.tcsh OUTPUT-DIR PARAM-SCRIPT BRANCH OBJECT(S)
##

if ($#argv != 4) then
  grep '^##' $0
  exit
endif

set outdir = $1
set params = $2
set branch = $3
set objects = ($4)

echo $branch
# read variables from input branch
source ./code/code.main/scripts-read-job-vars $branch "$objects" "genome genome_dir"

# run parameter script
source $params

# create path
scripts-create-path $outdir/

# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------

set hic_file = $branch/$objects/filtered.hic
set enzyme = `cut -f1-2,5-7 inputs/sample-sheet.tsv | grep -w "$objects" | cut -f4 | head -n1`

if ($tool == hint) then
  ./code/hicseq-translocations-hint.tcsh $outdir $params $genome $enzyme $hic_file
else
  echo "Error: Translocations tool $tool not supported." | scripts-send2err
endif

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------

# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."
