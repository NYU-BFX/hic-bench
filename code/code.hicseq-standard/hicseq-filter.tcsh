#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: hic-filter.tcsh OUTPUT-DIR PARAM-SCRIPT ALIGNMENT-BRANCH OBJECT(S)
##

# process command-line inputs
if ($#argv != 4) then
  grep '^##' $0
  exit
endif

set outdir = $1
set params = $2
set branch = $3
set objects = ($4)

# run parameter script
source $params

# indentify genome directory
set genome = `./code/read-sample-sheet.tcsh $sheet "$objects" genome | sort -u`
set enzyme = `./code/read-sample-sheet.tcsh $sheet "$objects" enzyme | sort -u`
set genome_dir = inputs/genomes/$genome


# create path
scripts-create-path $outdir/

if (-e $branch/$objects[1]/alignments.bam) then
#------------------------------------------------------------------------
# Case 1: a single alignments.bam file is available
#------------------------------------------------------------------------
  scripts-send2err "Filtering aligned reads..."
  
  set aligned_reads = `echo $objects | tr ' ' '\n' | awk -v d=$branch '{print d"/"$0"/alignments.bam"}'`
  if ($#objects == 1) then
    samtools view $aligned_reads | gtools-hic filter -v -E $genome_dir/$enzyme.fragments.bed --stats $outdir/stats_with_dups.tsv $filter_params | sort -t'	' -k2 >! $outdir/filtered_with_dups.reg
  else
    samtools merge - $aligned_reads | samtools view - | gtools-hic filter -v -E $genome_dir/$enzyme.fragments.bed --stats $outdir/stats_with_dups.tsv $filter_params | sort -t'	' -k2 >! $outdir/filtered_with_dups.reg
  endif

else if (-e $branch/$objects[1]/R1.bam) then
#------------------------------------------------------------------------
# Case 2: separate R1/R2 bam files are available
#------------------------------------------------------------------------
  echo -n '' >! $outdir/R12.sam
  foreach obj ($objects)
    scripts-send2err "Processing object $obj..."
    set objdir = $branch/$obj

    # process alignments
    scripts-send2err "-- Sorting R1 alignments..."
    samtools view $objdir/R1.bam | sort | tools-mergeuniq -merge -t '##%%%%%%%##' | sed 's/##%%%%%%%##.*//' >! $outdir/R1.sam
    scripts-send2err "-- Sorting R2 alignments..."
    samtools view $objdir/R2.bam | sort | tools-mergeuniq -merge -t '##%%%%%%%##' | sed 's/##%%%%%%%##.*//' >! $outdir/R2.sam
    scripts-send2err "-- Combining R1/R2 alignments into a single file..."
    sort -m $outdir/R1.sam $outdir/R2.sam >> $outdir/R12.sam

    # clean up
    scripts-send2err "-- Cleaning up..."
    rm -f $outdir/R[12].sam
  end
  
  # filtering
  scripts-send2err "Performing filtering using gtools..."
  cat $outdir/R12.sam | gtools-hic filter -v -E $genome_dir/$enzyme.fragments.bed --stats $outdir/stats_with_dups.tsv $filter_params | sort -t'	' -k2 >! $outdir/filtered_with_dups.reg
  rm -f $outdir/R12.sam

else if ($branch == 'inpdirs/align/results/align.by_sample.hicpro') then
#------------------------------------------------------------------------
# Case 3: his-pro allValidPairs files are available
#------------------------------------------------------------------------
  set allValidPairs = ()
  foreach obj ($objects)
    set allValidPairs = ($allValidPairs `ls -1 $branch/$obj/hic_results/data/$obj/*allValidPairs | grep -i allValidPairs`)
    #create a symbolic link to hic-pro output, used by loops whose inpdir is filter 
    ln -s ../../../../$branch/$obj/hic_results/data/$obj $outdir/hicpro
  end
  scripts-send2err "Converting hic-pro allValidPairs to filtered.reg format..."
  cat $allValidPairs | tools-cols -t 0 1 3 2 2 4 6 5 5 | tr '\t' ' ' | sed 's/ /\t/' >! $outdir/filtered.reg            # NOTE: can we filter by mapq here????
  set n_reads = `cat $outdir/filtered.reg | wc -l`
  set n_intra = `cat $outdir/filtered.reg | awk '$2==$6' | wc -l`
  set n_inter = $n_reads
  @ n_inter -= $n_intra
  set p_intra = `echo $n_intra/$n_reads | bc -l`
  set p_inter = `echo $n_inter/$n_reads | bc -l`
  gzip $outdir/filtered.reg
  ( echo "read-pairs $n_reads 1" ;\
    echo "unpaired 0 0" ;\
    echo "unmapped 0 0" ;\
    echo "multihit 0 0" ;\
    echo "single-sided 0 0" ;\
    echo "ds-no-fragment 0 0" ;\
    echo "ds-same-fragment 0 0" ;\
    echo "ds-too-close 0 0" ;\
    echo "ds-accepted-inter $n_inter $p_inter" ;\
    echo "ds-accepted-intra $n_intra $p_intra" ;\
    echo "ds-duplicate-inter 0 0" ;\
    echo "ds-duplicate-intra 0 0" ;\
    echo "ds-too-far 0 0" ;\
    echo "unclassified 0 0" ;\
  ) | tr ' ' '\t' >! $outdir/stats.tsv
  
  
else
#------------------------------------------------------------------------
# Case 4: validPairs files are available
#------------------------------------------------------------------------
  set valid_pairs = ()
  foreach obj ($objects)
    set valid_pairs = ($valid_pairs `ls -1 $branch/$obj/*gz | grep -i validpairs`)
  end
  scripts-send2err "Converting validPairs.txt to filtered.reg format..."
  cat $valid_pairs | gunzip | tools-cols -t 0 1 3 2 2 4 6 5 5 | tr '\t' ' ' | sed 's/ /\t/' >! $outdir/filtered.reg            # NOTE: can we filter by mapq here????
  set n_reads = `cat $outdir/filtered.reg | wc -l`
  set n_intra = `cat $outdir/filtered.reg | awk '$2==$6' | wc -l`
  set n_inter = $n_reads
  @ n_inter -= $n_intra
  set p_intra = `echo $n_intra/$n_reads | bc -l`
  set p_inter = `echo $n_inter/$n_reads | bc -l`
  gzip $outdir/filtered.reg
  ( echo "read-pairs $n_reads 1" ;\
    echo "unpaired 0 0" ;\
    echo "unmapped 0 0" ;\
    echo "multihit 0 0" ;\
    echo "single-sided 0 0" ;\
    echo "ds-no-fragment 0 0" ;\
    echo "ds-same-fragment 0 0" ;\
    echo "ds-too-close 0 0" ;\
    echo "ds-accepted-inter $n_inter $p_inter" ;\
    echo "ds-accepted-intra $n_intra $p_intra" ;\
    echo "ds-duplicate-inter 0 0" ;\
    echo "ds-duplicate-intra 0 0" ;\
    echo "ds-too-far 0 0" ;\
    echo "unclassified 0 0" ;\
  ) | tr ' ' '\t' >! $outdir/stats.tsv
  goto done
endif

if (-e $outdir/filtered_with_dups.reg) then 
  # remove duplicates
  scripts-send2err "Removing duplicates..."
  cat $outdir/filtered_with_dups.reg | uniq -f1 | gzip >! $outdir/filtered.reg.gz

  # update stats
  scripts-send2err "Updating statistics..."
  set n_intra_uniq = `cat $outdir/filtered_with_dups.reg | cut -f2 | uniq | cut -d' ' -f1,5 | awk '$1==$2' | wc -l`
  set n_inter_uniq = `cat $outdir/filtered_with_dups.reg | cut -f2 | uniq | cut -d' ' -f1,5 | awk '$1!=$2' | wc -l`
  Rscript ./code/update-filtered-stats.r $outdir/stats_with_dups.tsv $n_intra_uniq $n_inter_uniq >! $outdir/stats.tsv

  # cleanup
  rm -f $outdir/filtered_with_dups.reg $outdir/stats_with_dups.tsv
endif 

done:
# calculate distance statistics
./code/calc-distance-stats.tcsh $outdir/filtered.reg.gz >! $outdir/distance-stats.csv

# save variables
set >! $outdir/job.vars.tsv

scripts-send2err "Done."


