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
set field_names = (enzyme cell-type)              # field names to be included in the sample sheet header
set field_pos = (3 1)                             # field position in the sample name
set field_sep = (- -)                             # filed separator in the sample name
set n_fields = $#field_names

# create sample sheet
set inpdir = fastq
set sheet = sample-sheet.tsv
echo "sample group fastq-r1 fastq-r2 genome $field_names" | tr ' ' '\t' >! $sheet
foreach sample (`cd $inpdir; ls -1d *`)
  scripts-send2err "Importing sample $sample..."
  set group = `echo $sample | cut -d'-' -f-2`
  if (`ls -1 $inpdir/$sample | grep -c '.fastq.gz$'` == 0) then
    set fastq1 = '-'
	set fastq2 = '-'
  else
    set fastq1 = `cd $inpdir; ls -1 $sample/*_R1.fastq.gz $sample/*_R1_???.fastq.gz`
    set fastq2 = `cd $inpdir; ls -1 $sample/*_R2.fastq.gz $sample/*_R2_???.fastq.gz`
  endif
  set cell_type = `echo $sample | cut -d'-' -f1`
  echo -n "$sample\t$group\t`echo $fastq1 | tr ' ' ','`\t`echo $fastq2 | tr ' ' ','`\t$genome" >> $sheet
  set k = 1
  while ($k <= $n_fields) 
    set field = `echo $sample | cut -d$field_sep[$k] -f$field_pos[$k]`
    echo -n "\t$field" >> $sheet
	@ k ++
  end
  echo >> $sheet
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
echo "Configuring external data in data_external directory..."
echo "In this directory, you can include additional processed data, such are RNA-seq, ChIP-seq, ATAC-seq and other meta-data, organized by sample name, group name or cell type:"
mkdir -p data_external
foreach v (sample group cell-type)
#  echo "-- Creating directories for variable $v..."
  set k = `cat $sheet | head -1 | tr '\t' '\n' | grep -n "^$v"'$' | cut -d':' -f1`
  set X = `cat $sheet | cut -f$k | scripts-skipn 1 | sort -u`
  cd data_external
  mkdir -p $v
  cd $v
  mkdir -p $X
  cd ../..
end
echo




