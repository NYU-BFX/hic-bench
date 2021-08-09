#!/bin/tcsh
source ./code/code.main/custom-tcshrc         # customize shell environment

##
## USAGE: hicseq-align.tcsh OUTPUT-DIR PARAM-SCRIPT FASTQ-DIR OBJECT(S)
##

# process command-line inputs
if ($#argv != 4) then
  grep '^##' $0
  exit
endif

set out = $1
set params = $2
set fastq_dir = $3
set objects = ($4)

# test number of input objects
set object = $objects[1]
if ($#objects != 1) then
  scripts-send2err "Error: this operation allows only one input object!"
  exit 1
endif

# run parameter script
scripts-send2err "Setting parameters..."

source $params
if (! $?NSLOTS) then
  set threads = 16
else
  set threads = $NSLOTS
endif

# check if fastq files exist: if not, generate symbolic links to the input data files so that they can be used in the next step
if (`./code/read-sample-sheet.tcsh $sheet $object fastq-r1` == "-") then
  scripts-send2err "No input fastq files found, creating symbolic links to alternative input data."
  set object_dir = $fastq_dir/$object
  foreach f (`cd $object_dir; ls -1`)
   ln -s ../../../$object_dir/$f $out/$f
  end
  goto done
endif

# determine input fastq filenames
set fastq1 = `./code/read-sample-sheet.tcsh $sheet $object fastq-r1 | tr ',' '\n' | awk -v d=$fastq_dir '{print d"/"$0}'`
set fastq2 = `./code/read-sample-sheet.tcsh $sheet $object fastq-r2 | tr ',' '\n' | awk -v d=$fastq_dir '{print d"/"$0}'`

# align
if ($aligner == 'gtools') then                    ## Aligner = gtools

  if ( ($#fastq1 > 1) || ($#fastq2 > 1)) then
    scripts-send2err "Error: gtools-hic align does not allow multiple read pair files."
    exit
  endif
  gtools-hic align -v --work-dir $out/tmp -p $threads $align_params --reorder --bowtie-index $genome_index $fastq1 $fastq2 | samtools view -T $genome_index.fa -b1 - >! $out/alignments.bam
  rm -rf $out/tmp


else if ($aligner == 'bowtie2') then              ## Aligner = bowtie2

  set fastq1_comma = `echo $fastq1 | tr ' ' ','`
  set fastq2_comma = `echo $fastq2 | tr ' ' ','`
  set threads2 = `echo "$threads/1" | bc`

  # this code seems complicated, but it avoids generating intermediate sam files
  bash -c "paste -d'\n' <(bowtie2 -p $threads2 $align_params --reorder -x $genome_index --sam-nohead -U $fastq1_comma) <(bowtie2 -p $threads2 $align_params --reorder -x $genome_index --sam-nohead -U $fastq2_comma) | samtools view -T $genome_index.fa -b1 - > $out/alignments.bam"
  
  # old code below
  # bowtie2 -p $threads $align_params --reorder -x $genome_index --sam-nohead -U $fastq1_comma -S $out/alignments.R1.sam
  # bowtie2 -p $threads $align_params --reorder -x $genome_index --sam-nohead -U $fastq2_comma -S $out/alignments.R2.sam
  # paste -d'\n' $out/alignments.R1.sam $out/alignments.R2.sam | samtools view -T $genome_index.fa -b1 - >! $out/alignments.bam
  # rm -f $out/alignments.R1.sam $out/alignments.R2.sam

else if ($aligner == 'bwa') then              ## Aligner = bwa 

  # align R1 reads
  echo "Aligning R1 reads using BWA..." | scripts-send2err
  cat $fastq1 >! $out/R1.fastq.gz
  bwa mem $align_params -t $threads $genome_index $out/R1.fastq.gz | samtools view -Shb - >! $out/R1.bam
  rm -f $out/R1.fastq.gz
  
  # align R2 reads
  echo "Aligning R2 reads using BWA..." | scripts-send2err
  cat $fastq2 >! $out/R2.fastq.gz
  bwa mem $align_params -t $threads $genome_index $out/R2.fastq.gz | samtools view -Shb - >! $out/R2.bam
  rm -f $out/R2.fastq.gz
  
else if ($aligner == 'hicpro') then           ##run hic-pro 
  module unload r
  module load r/4.0.3
  module unload python
  module load python/cpu/3.7.2
  module load hic-pro

  set object_dir = $fastq_dir/$object

  #make custom config file from template 
  sed -e "s|refGenome|$genome|" \
      -e "s|enzyme_path|$PWD/$enzyme_path|" \
      -e "s|bowtie_idx_path|$PWD/$genome_index|" \
      -e "s|bowtie_param_global|$align_params_global|" \
      -e "s|bowtie_param_local|$align_params_local|" \
          $hicpro_config > config-hicpro-$objects.txt 
  
  
  
  mkdir $out/fastq/
  ln -s ../../../../$object_dir $out/fastq/
	
  #run hic-pro
  HiC-Pro -i $out/fastq/ -o $out -c config-hicpro-$objects.txt -s mapping -s proc_hic -s quality_checks -s merge_persample
  
  #clean up 
  rm -r $out/bowtie_results/bwt2_global $out/tmp
  rm $out/bowtie_results/bwt2/$objects/*.bwt2merged.bam
  

else
  scripts-send2err "Error: unknown aligner $aligner."
  exit
endif

# save variables
done:
set >! $out/job.vars.tsv




