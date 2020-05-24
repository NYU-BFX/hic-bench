#!/bin/tcsh

source ./inputs/params/params.tcsh

######## List of external files used for pygenometracks
# One can change the list of files by providing file paths that can be accessed.

# 1. hic_matrix
# The script automatically detects .h5 files from pipeline/tracks.
# 2. domains-diff
# The script automatically detects .bed file from pipeline/domains-diff.
# 3. compartments
# The script automatically detects .bedGraph files from pipeline/compartments.
# 4. virtual4C
# The script automatically detects .bedGraph.gz files from pipeline/virtual4C.

# E1. External .bigWig files (extension should end in .bw)
# The file link should be accessible from main-directory/pipelines/hicseq-standard/pipeline/pygenometracks-plot.
# Preferred order: ATAC-seq, H3K27ac ChIP-seq, RNA-seq
# Preferred color: #666666 for ATAC-seq, #FF33FF for H3K27ac ChIP-seq, #0033FF for RNA-seq

set external_bigwigs = ( \
  'inputs/data_external/group/ESC-untreated/ATAC-seq.ESC_treat_pileup.bw' \
  'inputs/data_external/group/MEF-untreated/ATAC-seq.MEF_treat_pileup.bw' \
  'inputs/data_external/group/ESC-untreated/H3K27ac-ChIP-seq.ESC_treat_pileup.bw' \
  'inputs/data_external/group/MEF-untreated/H3K27ac-ChIP-seq.MEF_treat_pileup.bw' \
  'inputs/data_external/group/ESC-untreated/RNA-seq.ESC.merged.scaled.bw' \
  'inputs/data_external/group/MEF-untreated/RNA-seq.MEF.merged.scaled.bw' \
)
set external_colors = ( \
  '#666666' \
  '#666666' \
  '#FF33FF' \
  '#FF33FF' \
  '#0033FF' \
  '#0033FF' \
)


# Final. Gene annotation
# The script automatically detects .gtf file from inputs/genome/$genome.

######## End of list

# Output settings
# A bed file may be set with the location accessible from pipelines/pygenometracks-plot.
# set inbed = params/Mtss1.and.Thada.bed
set inbed = "chr15:58,445,797-59,458,364"





# Now, process the directory structure of hic-bench for internal files.

set branch_medium = `echo $branch | sed 's/.*results\///'`
# $branch_medium example: tracks.by_group.h5.res_5kb/filter.by_sample.mapq_20/align.by_sample.bowtie2
# $branch_medium is what comes after inpdir/tracks/results
set group_var = `echo $branch_medium | cut -d'/' -f1 | cut -d'.' -f2`
# in running pygenometracks, $group_var is set to by_group, for now.

set branch_short = `echo $branch_medium | sed 's/.*filter\.by/filter.by/g'`
# $branch_short example: filter.by_sample.mapq_20/align.by_sample.bowtie2



# 1. hic_matrix
set library_size_branch = ../filter/results/$branch_short

# domains
set tad_caller = "hicratio.d_0500"
# tad_caller: "hicratio.d_0500" and "topdom.win5" are available

set domains_branch = ../domains/results/domains.$group_var.$tad_caller/matrix-ic.$group_var.cutoff_0/matrix-filtered.$group_var.res_40kb/$branch_short

# 2. domains-diff
# Select output from domains-diff
# is_normalize: "cpm" and "dist_norm" are available
# use_sample1_ref: TRUE or FALSE. TRUE selects "ref1" and FALSE selects "common".
set is_normalize = dist_norm
set use_sample1_ref = FALSE

if ( use_sample1_ref == TRUE ) then
  set tad_1_or_2 = ref1
else
  set tad_1_or_2 = common
endif

set domains_diff_branch = ../domains-diff/results/domains-diff.$group_var.$tad_caller.$is_normalize.$tad_1_or_2/matrix-ic.$group_var.cutoff_0/matrix-filtered.$group_var.res_40kb/$branch_short


# 3. compartments
set compartment_method = homer
set compartment_resolution = res_50kb

set compartment_dir = ../compartments/results/compartments.$group_var.$compartment_method.$compartment_resolution/$branch_short


# 4. virtual4C

set virtual4c_dir = ../virtual4C/results/virtual4C.by_group.win_10kb/matrix-sparse.by_group.unit_100bp.maxd_5Mb/$branch_short


### Finally, gene annotation goes here.





