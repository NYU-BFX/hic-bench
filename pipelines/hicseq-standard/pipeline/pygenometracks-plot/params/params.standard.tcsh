#!/bin/tcsh

source ./inputs/params/params.tcsh

set branch_medium = `echo $branch | sed 's/.*results\///'`
# $branch_medium example: tracks.by_group.h5.res_5kb/filter.by_sample.mapq_20/align.by_sample.bowtie2
# $branch_medium is what comes after inpdir/tracks/results
set group_var = `echo $branch_medium | cut -d'/' -f1 | cut -d'.' -f2`
# in running pygenometracks, $group_var is set to by_group, for now.

set branch_short = `echo $branch_medium | sed 's/.*filter\.by/filter.by/g'`
# $branch_short example: filter.by_sample.mapq_20/align.by_sample.bowtie2


# A bed file may be set with the exact location.
# set inbed = params/whgusdnflditlfh.bed
set inbed = ""
set outini = $outdir/hic_matrix.made.ini


# hic_matrix
set hic_matrix_file_1 = $inpdir/filtered.h5

set library_size_branch = ../filter/results/$branch_short
set scaling_factor_for_object = `Rscript ./code/hicseq-collect-library-size.r $library_size_branch $object`

# compartments
set compartment_method = homer
set compartment_resolution = res_50kb

set compartment_dir = ../compartments/results/compartments.$group_var.$compartment_method.$compartment_resolution/$branch_short

# domains
set tad_caller = "hicratio.d_0500"
# tad_caller: "hicratio.d_0500" and "topdom.win5" are available

set domains_branch = ../domains/results/domains.$group_var.$tad_caller/matrix-ic.$group_var.cutoff_0/matrix-filtered.$group_var.res_40kb/$branch_short

# domains-diff
set is_normalize = dist_norm
set use_sample1_ref = FALSE

if ( use_sample1_ref == TRUE ) then
  set tad_1_or_2 = ref1
else
  set tad_1_or_2 = common
endif

set domains_diff_branch = ../domains-diff/results/domains-diff.$group_var.$tad_caller.$is_normalize.$tad_1_or_2/matrix-ic.$group_var.cutoff_0/matrix-filtered.$group_var.res_40kb/$branch_short


# 1. hic_matrix
Rscript ./code/hicseq-prepare-ini-hicmatrix.r $hic_matrix_file_1 $object $scaling_factor_for_object $domains_branch/$object/domains.k=001.bed $outdir

# foreach line1 ( params/params.template.for.hic_matrix.txt )
#   if ( "$line1" == "file = template_and_modify_hic_matrix_file" ) then
#     echo $line1 | awk -v REP1="$hic_matrix_file_1" '{ gsub( /template_and_modify_hic_matrix_file/, REP1 ); print }' >> $outini
#   else if ( "$line1" == "title = template_and_modify_hic_matrix_title" ) then
#     echo $line1 | awk -v REP1="$object" '{ gsub( /template_and_modify_hic_matrix_title/, REP1 ); print }' >> $outini
#   else if ( "$line1" == "scale_factor = template_and_modify_hic_matrix_scale_factor" ) then
#     echo $line1 | awk -v REP1="$scaling_factor_for_object" '{ gsub( /template_and_modify_hic_matrix_scale_factor/, REP1 ); print }' >> $outini
#   else if ( "$line1" == "file = template_and_modify_domains_file" ) then
#     echo $line1 | awk -v REP1="$domains_branch\/$object\/domains.k=001.bed" '{ gsub( /template_and_modify_domains_file/, REP1 ); print }' >> $outini
#   else
#     cat $line1 >> $outini
#   endif

# end

# foreach line1 ( params/params.template.for.hic_matrix.txt )
#   if ( `echo "$line1" | grep -q "template_and_modify_hic_matrix_file"; echo $?` == 0 ) then
#     cat $line1 | sed "s/template_and_modify_hic_matrix_file/$hic_matrix_file_1/g" >> $outini
#   else if ( `echo "$line1" | grep -q "template_and_modify_hic_matrix_title"; echo $?` == 0 ) then
#     cat $line1 | sed "s/template_and_modify_hic_matrix_title/$object/g" >> $outini
#   else if ( `echo "$line1" | grep -q "template_and_modify_hic_matrix_scale_factor"; echo $?` == 0 ) then
#     cat $line1 | sed "s/template_and_modify_hic_matrix_scale_factor/$scaling_factor_for_object/g" >> $outini
#   else if ( `echo "$line1" | grep -q "template_and_modify_domains_file"; echo $?` == 0 ) then
#     cat $line1 | sed "s/template_and_modify_domains_file/$domains_branch\/$object\/domains.k=001.bed/g" >> $outini
#   else
#     cat $line1 >> $outini
#   endif

# end

# foreach line1 ( params/params.template.for.hic_matrix.txt )
#   if ( `echo "$line1" | grep -q "template_and_modify_hic_matrix_file"; echo $?` == 0 ) then
#     cat $line1 | sed "s/template_and_modify_hic_matrix_file/$hic_matrix_file_1/g" >> $outini
#   else if ( `echo "$line1" | grep -q "template_and_modify_hic_matrix_title"; echo $?` == 0 ) then
#     cat $line1 | sed "s/template_and_modify_hic_matrix_title/$object/g" >> $outini
#   else if ( `echo "$line1" | grep -q "template_and_modify_hic_matrix_scale_factor"; echo $?` == 0 ) then
#     cat $line1 | sed "s/template_and_modify_hic_matrix_scale_factor/$scaling_factor_for_object/g" >> $outini
#   else if ( `echo "$line1" | grep -q "template_and_modify_domains_file"; echo $?` == 0 ) then
#     cat $line1 | sed "s/template_and_modify_domains_file/$domains_branch\/$object\/domains.k=001.bed/g" >> $outini
#   else
#     cat $line1 >> $outini
#   endif

# end

# foreach line1 ( params/params.template.for.hic_matrix.txt )
#   if ( "$line1" =~ "*template_and_modify_hic_matrix_file" ) then
#     cat $line1 | sed "s/template_and_modify_hic_matrix_file/$hic_matrix_file_1/g" >> $outini
#   else if ( "$line1" =~ "*template_and_modify_hic_matrix_title" ) then
#     cat $line1 | sed "s/template_and_modify_hic_matrix_title/$object/g" >> $outini
#   else if ( "$line1" =~ "*template_and_modify_hic_matrix_scale_factor" ) then
#     cat $line1 | sed "s/template_and_modify_hic_matrix_scale_factor/$scaling_factor_for_object/g" >> $outini
#   else if ( "$line1" =~ "*template_and_modify_domains_file" ) then
#     cat $line1 | sed "s/template_and_modify_domains_file/$domains_branch\/$object\/domains.k=001.bed/g" >> $outini
#   else
#     cat $line1 >> $outini
#   endif

# end


# 2. domains-diff
# Select output from domains-diff
# is_normalize: "cpm" and "dist_norm" are available
# use_sample1_ref: TRUE or FALSE. TRUE selects "ref1" and FALSE selects "common".



# 3. compartments

# 4. virtual4C



### External data here.

# Preferred order of assays:

# ATAC-seq
# set external_atac = ./inputs/data_external/group/$object/ATAC-seq.ESC_treat_pileup.bw

# H3K27ac ChIP-seq
# set external_h3k27ac = ./inputs/data_external/group/$object/H3K27ac-ChIP-seq.ESC_treat_pileup.bw

# RNA-seq
# set external_rna = ./inputs/data_external/group/$object/RNA-seq.ESC.merged.scaled.bw

### Finally, gene annotation goes here.





