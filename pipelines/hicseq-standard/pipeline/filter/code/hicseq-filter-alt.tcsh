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
#----------------------------------------------
# Use gtools to filter
#----------------------------------------------
  scripts-send2err "Filtering aligned reads..."
  
  set aligned_reads = `echo $objects | tr ' ' '\n' | awk -v d=$branch '{print d"/"$0"/alignments.bam"}'`
  if ($#objects == 1) then
    samtools view $aligned_reads | gtools-hic filter -v -E $genome_dir/$enzyme.fragments.bed --stats $outdir/stats_with_dups.tsv $filter_params | sort -t'	' -k2 >! $outdir/filtered_with_dups.reg
  else
    samtools merge - $aligned_reads | samtools view - | gtools-hic filter -v -E $genome_dir/$enzyme.fragments.bed --stats $outdir/stats_with_dups.tsv $filter_params | sort -t'	' -k2 >! $outdir/filtered_with_dups.reg
  endif

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

else
#----------------------------------------------
# Use samtools to filter
#----------------------------------------------
  set N_reads = 0
  set N_mhits = 0
  set N_unmapped = 0
  set N_single = 0
  set N_close = 0
  set N_inter = 0
  set N_intra = 0
  set N_dups = 0
  echo -n '' >! $outdir/filtered.reg
  foreach obj ($objects)
    scripts-send2err "Processing object $obj..."
    set objdir = $branch/$obj

    # count number of read pairs
    scripts-send2err "-- Counting number of read pairs..."
    set n_reads = `samtools view $objdir/R1.bam | cut -f1 | sort -u | wc -l`

    # remove R1 multi-hits
    scripts-send2err "-- Removing multi-hits from R1 alignments..."	
    samtools view $objdir/R1.bam |  grep -v -e 'XA:Z:' -e 'SA:Z:' | gtools-regions reg | sort >! $outdir/R1_uniq.reg
    
    # remove R1 multi-hits
    scripts-send2err "-- Removing multi-hits from R2 alignments..."	
    samtools view $objdir/R2.bam |  grep -v -e 'XA:Z:' -e 'SA:Z:' | gtools-regions reg | sort >! $outdir/R2_uniq.reg
    
    # merge R1 and R2
    scripts-send2err "-- Merging unique R1/R2 alignments..."
    sort -m $outdir/R1_uniq.reg $outdir/R2_uniq.reg | tools-mergeuniq -merge >! $outdir/filtered_uniq.reg
    set nf_uniq = `cat $outdir/filtered_uniq.reg | wc -l`
    set n_mhits = `echo $n_reads-$nf_uniq | bc`
    #rm -f $outdir/R*_uniq.reg
    
    # remove unmapped
    scripts-send2err "-- Removing unmapped read pairs..."
    cat $outdir/filtered_uniq.reg | awk '($2!="*") && ($4!="*")' >! $outdir/filtered_unmapped.reg
    set nf_unmapped = `cat $outdir/filtered_unmapped.reg | wc -l`
    set n_unmapped = `echo $nf_uniq-$nf_unmapped | bc`
    
    # remove single-sided
    scripts-send2err "-- Removing single-sided alignments..."
    cat $outdir/filtered_unmapped.reg | awk '$6!=""' >! $outdir/filtered_2sided.reg
    set nf_2sided = `cat $outdir/filtered_2sided.reg | wc -l`
    set n_single = `echo $nf_unmapped-$nf_2sided | bc`

    # remove PCR duplicates
    scripts-send2err "-- Removing PCR duplicates..."
    cat $outdir/filtered_2sided.reg | sort -t'	' -k2 | uniq -f1 >! $outdir/filtered_nodups.reg 
    set nf_nodups = `cat $outdir/filtered_nodups.reg | wc -l `
    set n_dups = `echo $nf_2sided-$nf_nodups | bc`

    # remove reads pairs that are too close (FINAL)
    scripts-send2err "-- Removing reads pairs that are too close..."
    cat $outdir/filtered_nodups.reg | awk -v d=$mindist '($2!=$6) || ($4-$8>=d) || ($8-$4>=d)' >! $outdir/filtered_notclose.reg            
    set nf_notclose = `cat $outdir/filtered_notclose.reg | wc -l `
    set n_close = `echo $nf_nodups-$nf_notclose | bc`

    # count intra/inter read pairs
    scripts-send2err "-- Counting usable intra/inter-chromosomal read pairs..."
    set n_inter = `cat $outdir/filtered_notclose.reg | awk '$2!=$6' | wc -l`
    set n_intra = `cat $outdir/filtered_notclose.reg | awk '$2==$6' | wc -l`

    # add final filtered reads to main file
    scripts-send2err "-- Adding filtered reads to filtered.reg file..."
    cat $outdir/filtered_notclose.reg >> $outdir/filtered.reg
	
    # clean up
    scripts-send2err "-- Cleaning up..."
    #rm -f $outdir/filtered_*.reg

    # update stats
    @ N_reads += $n_reads
    @ N_mhits += $n_mhits
    @ N_unmapped += $n_unmapped
    @ N_single += $n_single
    @ N_close += $n_close
    @ N_inter += $n_inter
    @ N_intra += $n_intra
    @ N_dups += $n_dups
  end
  
  # gzip filtered reads
  scripts-send2err "Gzipping filtered.reg file..."
  gzip $outdir/filtered.reg

  # store stats
  (echo "read-pairs	$N_reads	1"; \
   echo "unpaired	0	0" ; \
   echo "unmapped	$N_unmapped	`echo $N_unmapped/$N_reads | bc -l`" ; \
   echo "multihit	$N_mhits	`echo $N_mhits/$N_reads | bc -l`" ; \
   echo "single-sided	$N_single	`echo $N_single/$N_reads | bc -l`" ; \
   echo "ds-no-fragment	0	0" ; \
   echo "ds-same-fragment	0	0" ; \
   echo "ds-too-close	$N_close	`echo $N_close/$N_reads | bc -l`" ; \
   echo "ds-accepted-inter	$N_inter	`echo $N_inter/$N_reads | bc -l`" ; \
   echo "ds-accepted-intra	$N_intra	`echo $N_intra/$N_reads | bc -l`" ; \
   echo "ds-duplicate-inter	0	0" ; \
   echo "ds-duplicate-intra	$N_dups	`echo $N_dups/$N_reads | bc -l`" ; \
   echo "ds-too-far	0	0" ; \
   echo "unclassified	0	0" ; \
  ) >! $outdir/stats.tsv

endif


# save variables
set >! $outdir/job.vars.tsv

scripts-send2err "Done."


