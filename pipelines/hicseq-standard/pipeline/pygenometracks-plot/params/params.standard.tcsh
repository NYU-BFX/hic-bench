#!/bin/tcsh

source ./inputs/params/params.tcsh

### Instructions!
# Wherever USER INPUT is seen, the code is meant to be replaced to reflect the characteristics of directory structure.



# First, process the directory structure of hic-bench for internal files.

set branch_medium = `echo $branch | sed 's/.*results\///'`
# $branch_medium example: tracks.by_group.h5.res_5kb/filter.by_sample.mapq_20/align.by_sample.bowtie2
# $branch_medium is what comes after inpdir/tracks/results
set group_var = `echo $branch_medium | cut -d'/' -f1 | cut -d'.' -f2`
# in running pygenometracks, $group_var is set to by_group, for now.

set branch_short = `echo $branch_medium | sed 's/.*filter\.by/filter.by/g'`
# $branch_short example: filter.by_sample.mapq_20/align.by_sample.bowtie2


# Information of input files is written below.

# 1. hic_matrix
# Select directory
# $branch_short (pipeline/tracks) is selected as the hic_matrix data directory by default.
# The script automatically detects .h5 files from $branch_short.

# Access library size
set library_size_branch = ../filter/results/$branch_short

# Access domains calls
### USER INPUT!
set tad_caller = "hicratio.d_0500"
# tad_caller: "hicratio.d_0500" and "topdom.win5" are available
### USER INPUT! end

set domains_branch = ../domains/results/domains.$group_var.$tad_caller/matrix-ic.$group_var.cutoff_0/matrix-filtered.$group_var.res_40kb/$branch_short

# 2. domains-diff
# Select output from domains-diff

### USER INPUT!
set is_normalize = dist_norm
set use_sample1_ref = FALSE
# is_normalize: "cpm" and "dist_norm" are available
# use_sample1_ref: TRUE or FALSE. TRUE selects "ref1" and FALSE selects "common".
### USER INPUT! end

if ( use_sample1_ref == TRUE ) then
  set tad_1_or_2 = ref1
else
  set tad_1_or_2 = common
endif

set domains_diff_branch = ../domains-diff/results/domains-diff.$group_var.$tad_caller.$is_normalize.$tad_1_or_2/matrix-ic.$group_var.cutoff_0/matrix-filtered.$group_var.res_40kb/$branch_short

### USER INPUT!
# Choose $object_pair among the directory name inside $domains_diff_branch.
set object_pair = "ESC-untreated.MEF-untreated"
### USER INPUT! end


# Finally, the script detects .bed file from pipeline/domains-diff.

# 3. compartments
# Select directory

### USER INPUT!
set compartment_method = homer
set compartment_resolution = res_50kb
### USER INPUT! end

set compartment_dir = ../compartments/results/compartments.$group_var.$compartment_method.$compartment_resolution/$branch_short

# The script automatically detects .bedGraph files from pipeline/compartments.

# 4. virtual4C
# Select directory
set virtual4C_dir = ../virtual4C/results/virtual4C.by_group.win_10kb/matrix-sparse.by_group.unit_100bp.maxd_5Mb/$branch_short

### USER INPUT!
# One must select the exact .bedGraph files for each run for a region.
# Usual instructions for using single quotes versus double quotes (string interpolation) apply.
set v4c_bdgs = ( \
  "$virtual4C_dir/ESC-untreated/Mtss1-chr15-v4C.bedgraph.gz" \
  "$virtual4C_dir/MEF-untreated/Mtss1-chr15-v4C.bedgraph.gz" \
)
set v4c_colors = ( \
  'darkred' \
  'darkblue' \
)
### USER INPUT!

# E1. External .bigWig files (extension should end in .bw)
# The file link should be accessible from main-directory/pipelines/hicseq-standard/pipeline/pygenometracks-plot.
# Preferred order: ATAC-seq, H3K27ac ChIP-seq, RNA-seq
# Preferred color: #666666 for ATAC-seq, #FF33FF for H3K27ac ChIP-seq, #0033FF for RNA-seq

### USER INPUT!
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
### USER INPUT! end

# Final. Gene annotation
# The script automatically detects .gtf file from inputs/genome/$genome.
set GENE_GTF = '/gpfs/data/abl/home/choh09/iGenomes/Mus_musculus/UCSC/mm10/Annotation/Genes/gencode.vM25.annotation.gtf'
# This part not complete
# set GENE_GTF = "$genome_dir/gencode.vM25.annotation.gtf"

######## End of list

# Output settings
# A bed file may be set with the location accessible from pipelines/pygenometracks-plot.
# set inbed = params/Mtss1.and.Thada.bed
set inbed = "chr15:58,445,797-59,458,364"










### Finally, gene annotation goes here.





