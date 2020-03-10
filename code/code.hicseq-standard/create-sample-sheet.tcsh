#!/bin/tcsh
source ./code/code.main/custom-tcshrc    # customize shell environment

##
## USAGE: create-sample-sheet.tcsh GENOME={hg19,mm10}
##
## FUNCTION: create sample sheet automatically from input files in fastq directory
##

# process command-line inputs
if ($#argv != 1) then
  grep '^##' $0
  exit
endif

set genome = $1
set enzyme_field = 3

# create sample sheet
set inpdir = fastq
set sheet = sample-sheet.tsv
echo "sample group fastq-r1 fastq-r2 genome enzyme cell-type" | tr ' ' '\t' >! $sheet
foreach sample (`cd $inpdir; ls -1d *`)
  scripts-send2err "Importing sample $sample..."
  set group = `echo $sample | cut -d'-' -f-2`
  set fastq1 = `cd $inpdir; ls -1 $sample/*_R1.fastq.gz $sample/*_R1_???.fastq.gz`
  set fastq2 = `cd $inpdir; ls -1 $sample/*_R2.fastq.gz $sample/*_R2_???.fastq.gz`
  set enzyme = `echo $sample | cut -d'-' -f$enzyme_field`
  set cell_type = `echo $sample | cut -d'-' -f1`
  echo "$sample\t$group\t`echo $fastq1 | tr ' ' ','`\t`echo $fastq2 | tr ' ' ','`\t$genome\t$enzyme\t$cell_type" >> $sheet
end

echo "Your sample sheet has been created! Here is how it looks:"
echo
cat $sheet
echo

echo "Diagnostics: "
echo "Field #1: "
cat $sheet | scripts-skipn 1 | cut -f1 | cut -d'-' -f1 | sort | uniq -c
echo "Field #2: "
cat $sheet | scripts-skipn 1 | cut -f1 | cut -d'-' -f2 | sort | uniq -c
echo "Field #3: "
cat $sheet | scripts-skipn 1 | cut -f1 | cut -d'-' -f3 | sort | uniq -c

echo 
echo "You can configure processed external data (e.g. ChIP-seq, ATAC-seq) in the data_external directory. Please read the README file. "
echo
