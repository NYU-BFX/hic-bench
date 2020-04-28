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
set branch = $3
set object1 = $4
set object2 = $5

set objects = ($object1 $object2)
set inpdir1 = $branch/$object1

# read variables from input branch
source ./code/code.main/scripts-read-job-vars $branch "$objects" "genome genome_dir bin_size"

# run parameter script
source $params

# create path
scripts-create-path $outdir/

# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------


if ($tool == "tadsplimer") then
  scripts-send2err "Running TADsplimer for every chromosome."
  # run TADsplimer for every chromosome
  set jid = 

  # This code assumes $inpdir1 and $inpdir2 have the same chromosome structure. Be careful
  set matrices = `cd $inpdir1; ls -1 matrix.*.tsv | grep -vwE "$chrom_excluded"`
  foreach mat ($matrices)
    scripts-send2err "Processing $mat..."

    # Actual run script
    set jpref = $outdir/__jdata/job.`echo $mat | sed 's/\.[^.]\+$//'`
    scripts-create-path $jpref
    set jid = ($jid `scripts-qsub-run $jpref 1 "16G" ./code/insulation-diff-tadsplimer.sh $outdir $branch $mat $object1 $object2`)
  end

endif

# wait until all jobs are completed
scripts-send2err "Waiting until all jobs are completed..."
scripts-qsub-wait "$jid"


# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------

# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."


