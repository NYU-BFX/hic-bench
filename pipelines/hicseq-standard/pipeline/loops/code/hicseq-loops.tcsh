#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: hicseq-compartments.tcsh OUTPUT-DIR PARAM-SCRIPT BRANCH OBJECT(S)
##

if ($#argv != 4) then
  grep '^##' $0
  exit
endif

set outdir = $1
set params = $2
set branch = $3
set objects = ($4)

# read variables from input branch
source ./code/code.main/scripts-read-job-vars $branch "$objects" "genome genome_dir"

# run parameter script
source $params

# create path
scripts-create-path $outdir/

# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------

# make a list of filtered read files for all input objects
set reg_files = `echo $objects | tr ' ' '\n' | awk -v d=$branch '{print d"/"$0"/filtered.reg.gz"}'`

if ($tool == fithic) then
	./code/hicseq-loops-fithic.tcsh $outdir $params "$reg_files" $genome
else
  	echo "Error: Compartment calling tool $tool not supported." | scripts-send2err
endif


# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------

# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."
