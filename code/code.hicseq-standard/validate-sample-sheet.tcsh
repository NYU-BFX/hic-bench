#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # custom shell environment

##
## USAGE: validate-sample-sheet SAMPLE-SHEET-TSV
##

if ($#argv != 1) then
  grep '^##' $0
  exit 1
endif

set sheet = $1

# check for SPACES
if (`cat $sheet | grep ' ' | wc -l` > 0) then
  scripts-send2err 'Error: spaces are not allowed in sample sheet.'
  exit 1
endif

# duplicate sample names (column #1)
if (`cat $sheet | scripts-skipn 1 | cut -f1 | sort | uniq -d | wc -l` > 0) then
  scripts-send2err 'Error: duplicate sample names are not allowed in column 1 of the sample sheet.'
  exit 1
endif

# file names should exist in fastq directory
set files = `cat $sheet | scripts-skipn 1 | cut -f3,4 | tr ',\t' '\n' | sort -u`
foreach f ($files)
  if (! -e inputs/fastq/$f && $f != -) then
    scripts-send2err "Error: file \'$f\' not found in inputs/fastq directory."
    exit 1
  endif
end



