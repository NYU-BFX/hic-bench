#!/bin/tcsh
source ./code/code.main/custom-tcshrc    # customize shell environment

##
## USAGE: hicseq-matrix-diff.tcsh OUTPUT-DIR PARAM-SCRIPT MATRIX-BRANCH OBJECT1 OBJECT2
##

if ($#argv != 5) then
ste(obj1, "specific SE", sep = " ")
  grep '^##' $0
  exit
endif

set outdir = $1
set params = $2
set branch = $3
set object1 = $4
set object2 = $5

set objects = ($object1 $object2)

# read variables from input branch
source ./code/code.main/scripts-read-job-vars $branch "$objects" "genome genome_dir bin_size"

# run parameter script
source $params

# create path
scripts-create-path $outdir/

# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------

# code dir
set codedir = code/hicseq-domains-diff-scripts

# First, assign TAD call files for the samples of interest
set domains1 = $domains_branch/$object1/domains.k=001.bed
set domains2 = $domains_branch/$object2/domains.k=001.bed

if ($perform_analysis == TRUE) then
	if ($use_sample1_ref == TRUE) then
		# Use sample 1 TADs as reference
		# Perform Hi-C fold-change analysis
		cat $domains1 >! $outdir/domains1.tsv
		cat $domains1 >! $outdir/domains2.tsv
		cat $domains1 >! $outdir/domains_common.tsv
	else
		# Find consistent TADs between samples
		perl $codedir/find-consistent-domains.pl $domains1 $domains2 $max_boundary_dist $bin_size $outdir/domains1.tsv $outdir/domains2.tsv $outdir/domains_common.tsv
	endif

	# Perform Hi-C fold-change analysis
	$codedir/run_comparison.sh $codedir/differential_tad_activity.r $branch/$object2 $branch/$object1 $outdir/domains2.tsv $outdir/domains1.tsv $min_tad_size $max_tad_size $is_normalize $centrotelo_file $bin_size $outdir

	# Create boxplot with differentially active TADs
	R --no-save $outdir/final_results.tsv $gene_name FALSE FALSE \
		$object2 $object1 $bin_size $min_tad_size $outdir/final_results < $codedir/differential_tad_activity_expression.r
endif

# Integrate RNA-Seq data
if ($rnaseq == TRUE) then
 
	set rna_branch = ../domains-diff-integration/rnaseq
	set rna1 = $rna_branch/$object1/rnaseq_cpm.rds
	set rna2 = $rna_branch/$object2/rnaseq_cpm.rds

	# Add boxplot plot with the RNA-Seq integration
	R --no-save $outdir/final_results.tsv $gene_name $rna1 $rna2 \
        	$object2 $object1 $bin_size $min_tad_size $outdir/final_results < $codedir/differential_tad_activity_expression.r
endif

# Integrate superenhancer data
# Superenhancer calling must have been performed prior to this analysis 
# and the files need to be saved in the corresponding directory as
# indicated by the following commands.
if ($superenhancers == TRUE) then
	set analysis_type = superenhancer
	set se_branch = ../domains-diff-integration/superenhancers
	set se1 = $se_branch/$object1/superenhancers.bed
	set se2 = $se_branch/$object2/superenhancers.bed
	
	# Overlap with all common TADs
	bedtools intersect -a $outdir/domains_common.tsv -b $se1 -wa -wb > $outdir/all_TADs_overlap_se1.tsv
	bedtools intersect -a $outdir/domains_common.tsv -b $se2 -wa -wb > $outdir/all_TADs_overlap_se2.tsv
	
	# Count the number of common TADs and SE peaks overlapping them
	set num_common_TADs=`wc -l $outdir/domains_common.tsv | awk '{print $1}'`	
	
	set num_peaks1=`wc -l $outdir/all_TADs_overlap_se1.tsv | awk '{print $1}'`	
	set num_peaks2=`wc -l $outdir/all_TADs_overlap_se2.tsv | awk '{print $1}'`	
	
	# Note. Active TADs are the ones more active in object1
	bedtools intersect -a $outdir/final_results_active-TADs.bed -b $se1 -wa -wb > $outdir/active_TADs_overlap_se1.tsv
	bedtools intersect -a $outdir/final_results_active-TADs.bed -b $se2 -wa -wb > $outdir/active_TADs_overlap_se2.tsv

	bedtools intersect -a $outdir/final_results_inactive-TADs.bed -b $se1 -wa -wb > $outdir/inactive_TADs_overlap_se1.tsv
	bedtools intersect -a $outdir/final_results_inactive-TADs.bed -b $se2 -wa -wb > $outdir/inactive_TADs_overlap_se2.tsv

	bedtools intersect -a $outdir/final_results_unchanged-TADs.bed -b $se1 -wa -wb > $outdir/unchanged_TADs_overlap_se1.tsv
	bedtools intersect -a $outdir/final_results_unchanged-TADs.bed -b $se2 -wa -wb > $outdir/unchanged_TADs_overlap_se2.tsv
	
	# Plot the number of Superenhancers per TAD for TADs in different categories
	# This is the same script used for Approach #1 and #2 of ATAC-Seq data integration
	R --no-save 	$outdir/active_TADs_overlap_se1.tsv $outdir/active_TADs_overlap_se2.tsv \
			$outdir/inactive_TADs_overlap_se1.tsv $outdir/inactive_TADs_overlap_se2.tsv \
			$outdir/unchanged_TADs_overlap_se1.tsv $outdir/unchanged_TADs_overlap_se2.tsv \
                	$object2 $object1 $bin_size $outdir/final_results $analysis_type \
			$num_peaks1 $num_peaks2 $num_common_TADs $outdir/domains_common.tsv \
			$outdir/final_results_active-TADs.bed $outdir/final_results_inactive-TADs.bed \
			$outdir/final_results_unchanged-TADs.bed  < $codedir/differential_number_ChIP_Seq_peaks.r
	
endif

# Integrate typical enhancer data
if ($enhancers == TRUE) then
	set analysis_type = enhancers
	set enhancer_branch = ../domains-diff-integration/enhancers
	set enhancers1 = $enhancer_branch/$object1/H3K27ac_peaks.bed
	set enhancers2 = $enhancer_branch/$object2/H3K27ac_peaks.bed
	
	# Remove promoter peaks - those overlapping TSSs
	perl $codedir/split_peaks_coding-vs-noncoding.pl $enhancers1 $genome_file $outdir/promoter_peaks_enhancers1.bed $outdir/non_promoter_peaks_enhancers1.bed
	perl $codedir/split_peaks_coding-vs-noncoding.pl $enhancers2 $genome_file $outdir/promoter_peaks_enhancers2.bed $outdir/non_promoter_peaks_enhancers2.bed	

	# Integrate with differentially active TADs
	bedtools intersect -a $outdir/domains_common.tsv -b $outdir/non_promoter_peaks_enhancers1.bed -wa -wb > $outdir/all_TADs_overlap_enhancers1.tsv
	bedtools intersect -a $outdir/domains_common.tsv -b $outdir/non_promoter_peaks_enhancers2.bed -wa -wb > $outdir/all_TADs_overlap_enhancers2.tsv
	
	# Count the number of common TADs
	set num_common_TADs=`wc -l $outdir/domains_common.tsv | awk '{print $1}'`	
	
	# Count the enhancer peaks overlapping common TADs
	set num_peaks1=`wc -l $outdir/all_TADs_overlap_enhancers1.tsv | awk '{print $1}'`	
	set num_peaks2=`wc -l $outdir/all_TADs_overlap_enhancers2.tsv | awk '{print $1}'`	
	
	bedtools intersect -a $outdir/final_results_active-TADs.bed -b $outdir/non_promoter_peaks_enhancers1.bed -wa -wb > $outdir/active_TADs_overlap_enhancers1.tsv
	bedtools intersect -a $outdir/final_results_active-TADs.bed -b $outdir/non_promoter_peaks_enhancers2.bed -wa -wb > $outdir/active_TADs_overlap_enhancers2.tsv

	bedtools intersect -a $outdir/final_results_inactive-TADs.bed -b $outdir/non_promoter_peaks_enhancers1.bed -wa -wb > $outdir/inactive_TADs_overlap_enhancers1.tsv
	bedtools intersect -a $outdir/final_results_inactive-TADs.bed -b $outdir/non_promoter_peaks_enhancers2.bed -wa -wb > $outdir/inactive_TADs_overlap_enhancers2.tsv

	bedtools intersect -a $outdir/final_results_unchanged-TADs.bed -b $outdir/non_promoter_peaks_enhancers1.bed -wa -wb > $outdir/unchanged_TADs_overlap_enhancers1.tsv
	bedtools intersect -a $outdir/final_results_unchanged-TADs.bed -b $outdir/non_promoter_peaks_enhancers2.bed -wa -wb > $outdir/unchanged_TADs_overlap_enhancers2.tsv

	# Plot the number of Superenhancers per TAD for TADs in different categories
	# This is the same script used for Approach #1 and #2 of ATAC-Seq data integration
	R --no-save 	$outdir/active_TADs_overlap_enhancers1.tsv $outdir/active_TADs_overlap_enhancers2.tsv \
			$outdir/inactive_TADs_overlap_enhancers1.tsv $outdir/inactive_TADs_overlap_enhancers2.tsv \
			$outdir/unchanged_TADs_overlap_enhancers1.tsv $outdir/unchanged_TADs_overlap_enhancers2.tsv \
                	$object2 $object1 $bin_size $outdir/final_results $analysis_type \
			$num_peaks1 $num_peaks2 $num_common_TADs $outdir/domains_common.tsv \
			$outdir/final_results_active-TADs.bed $outdir/final_results_inactive-TADs.bed \
			$outdir/final_results_unchanged-TADs.bed  < $codedir/differential_number_ChIP_Seq_peaks.r

endif

# Integrate any kind of ChIP data - ATAC Seq in this senario
if ($atacseq == TRUE) then
	set analysis_type = ATAC_Seq
	set chip_branch = ../domains-diff-integration/ATAC-Seq
	set diff_bind_file = $chip_branch/diff_bind."$object2"-vs-"$object1".q100.csv
	set chip1 = $chip_branch/$object1/peaks.bed
	set chip2 = $chip_branch/$object2/peaks.bed

	# APPROACH #1
	# Integrate total number of ChIP-Seq peaks with differentially active TADs
	
	bedtools intersect -a $outdir/domains_common.tsv -b $chip1 -wa -wb > $outdir/all_TADs_overlap_chip1.tsv
	bedtools intersect -a $outdir/domains_common.tsv -b $chip2 -wa -wb > $outdir/all_TADs_overlap_chip2.tsv
	
	# Count the number of common TADs
	set num_common_TADs=`wc -l $outdir/domains_common.tsv | awk '{print $1}'`	
	
	# Count total number of peaks per object for normalization
	set num_peaks1=`wc -l $outdir/all_TADs_overlap_chip1.tsv | awk '{print $1}'`	
	set num_peaks2=`wc -l $outdir/all_TADs_overlap_chip2.tsv | awk '{print $1}'`	
	
	# Count number of peaks for each differentially active TAD across the two samples
	bedtools intersect -a $outdir/final_results_active-TADs.bed -b $chip1 -wa -wb > $outdir/active_TADs_overlap_chip1.tsv
        bedtools intersect -a $outdir/final_results_active-TADs.bed -b $chip2 -wa -wb > $outdir/active_TADs_overlap_chip2.tsv

        bedtools intersect -a $outdir/final_results_inactive-TADs.bed -b $chip1 -wa -wb > $outdir/inactive_TADs_overlap_chip1.tsv
        bedtools intersect -a $outdir/final_results_inactive-TADs.bed -b $chip2 -wa -wb > $outdir/inactive_TADs_overlap_chip2.tsv

        bedtools intersect -a $outdir/final_results_unchanged-TADs.bed -b $chip1 -wa -wb > $outdir/unchanged_TADs_overlap_chip1.tsv
        bedtools intersect -a $outdir/final_results_unchanged-TADs.bed -b $chip2 -wa -wb > $outdir/unchanged_TADs_overlap_chip2.tsv

	# Plot the number of ATAC-Seq peaks per TAD for TADs in different categories
	R --no-save 	$outdir/active_TADs_overlap_chip1.tsv $outdir/active_TADs_overlap_chip2.tsv \
			$outdir/inactive_TADs_overlap_chip1.tsv $outdir/inactive_TADs_overlap_chip2.tsv \
			$outdir/unchanged_TADs_overlap_chip1.tsv $outdir/unchanged_TADs_overlap_chip2.tsv \
                	$object2 $object1 $bin_size $outdir/final_results $analysis_type  \
			$num_peaks1 $num_peaks2 $num_common_TADs $outdir/domains_common.tsv \
			$outdir/final_results_active-TADs.bed $outdir/final_results_inactive-TADs.bed \
			$outdir/final_results_unchanged-TADs.bed  < $codedir/differential_number_ChIP_Seq_peaks.r

	# APPROACH #2
	# Integrate the differential peaks - the ones which are stronger to each condition	
	set analysis_type = diff_ATAC_Seq
	
	# Determine peaks stronger in obj1 and obj2
	# I am assuming that the fold calculation will be calculated as obj1 - obj2 
	# Filtering for significant changes: FDR 0.05
	awk -F"," '{ if (($9 < 0)&&($11 < 0.05)) { print } }' $diff_bind_file | sed 's/,/\t/g' | sed 's/"//g' | tail -n +2  > $outdir/diff_bind_higher_"$object1".csv 
	awk -F"," '{ if (($9 > 0)&&($11 < 0.05)) { print }}' $diff_bind_file  | sed 's/,/\t/g' | sed 's/"//g' | tail -n +2  > $outdir/diff_bind_higher_"$object2".csv 
	
	bedtools intersect -a $outdir/domains_common.tsv -b $outdir/diff_bind_higher_"$object2".csv -wa -wb > $outdir/all_TADs_overlap_chip1.tsv
	bedtools intersect -a $outdir/domains_common.tsv -b $outdir/diff_bind_higher_"$object1".csv -wa -wb > $outdir/all_TADs_overlap_chip2.tsv
	
	# Count the number of common TADs
	set num_common_TADs=`wc -l $outdir/domains_common.tsv | awk '{print $1}'`	
	
	# Count total number of peaks per object for normalization
	set num_peaks1=`wc -l $outdir/all_TADs_overlap_chip1.tsv | awk '{print $1}'`	
	set num_peaks2=`wc -l $outdir/all_TADs_overlap_chip2.tsv | awk '{print $1}'`	

	# Overlap the peaks with the differentially active TADs
	bedtools intersect -a $outdir/final_results_active-TADs.bed -b $outdir/diff_bind_higher_"$object2".csv -wa -wb > $outdir/active_TADs_overlap_diff_bind_chip1.tsv
        bedtools intersect -a $outdir/final_results_active-TADs.bed -b $outdir/diff_bind_higher_"$object1".csv -wa -wb > $outdir/active_TADs_overlap_diff_bind_chip2.tsv

        bedtools intersect -a $outdir/final_results_inactive-TADs.bed -b $outdir/diff_bind_higher_"$object2".csv -wa -wb > $outdir/inactive_TADs_overlap_diff_bind_chip1.tsv
        bedtools intersect -a $outdir/final_results_inactive-TADs.bed -b $outdir/diff_bind_higher_"$object1".csv -wa -wb > $outdir/inactive_TADs_overlap_diff_bind_chip2.tsv

        bedtools intersect -a $outdir/final_results_unchanged-TADs.bed -b $outdir/diff_bind_higher_"$object2".csv -wa -wb > $outdir/unchanged_TADs_overlap_diff_bind_chip1.tsv
        bedtools intersect -a $outdir/final_results_unchanged-TADs.bed -b $outdir/diff_bind_higher_"$object1".csv -wa -wb > $outdir/unchanged_TADs_overlap_diff_bind_chip2.tsv

	# Plot the number of differential ATAC-Seq peaks per TAD for TADs in different categories
	R --no-save 	$outdir/active_TADs_overlap_diff_bind_chip1.tsv $outdir/active_TADs_overlap_diff_bind_chip2.tsv \
			$outdir/inactive_TADs_overlap_chip1.tsv $outdir/inactive_TADs_overlap_chip2.tsv \
			$outdir/unchanged_TADs_overlap_chip1.tsv $outdir/unchanged_TADs_overlap_chip2.tsv \
                	$object2 $object1 $bin_size $outdir/final_results $analysis_type \
			$num_peaks1 $num_peaks2 $num_common_TADs $outdir/domains_common.tsv \
			$outdir/final_results_active-TADs.bed $outdir/final_results_inactive-TADs.bed \
			$outdir/final_results_unchanged-TADs.bed  < $codedir/differential_number_ChIP_Seq_peaks.r

	# APPROACH #3
	# Calculate the logFC of ATAC-Seq peaks that are differentially strong across the two samples
	# as identified by ChIP-Seq analysis and diff_bind() and plot the logFC distribution.
	
	# For this task we need to extend the TAD calls based on the TAD caller to account for mistakes
	# Extension should be 0 for crane TAD calling and 80kb for hicratio
	#
	set extension = 0
	if ($tad_caller == hicratio.d_0500) then 
		set extension = 80000
	endif

	echo $tad_caller
	echo $extension

	perl $codedir/overlap_intra_TAD_ChIP_Seq_peaks.pl \
        	$outdir/final_results_active-TADs.tsv \
        	$outdir/diff_bind."$object2"-vs-"$object1".q100-zscore.bed \
        	$extension \
        	$outdir/active-vs-diffBind-FC_intra_TAD.tsv 

	perl $codedir/overlap_intra_TAD_ChIP_Seq_peaks.pl \
        	$outdir/final_results_inactive-TADs.tsv \
        	$outdir/diff_bind."$object2"-vs-"$object1".q100-zscore.bed \
        	$extension \
        	$outdir/inactive-vs-diffBind-FC_intra_TAD.tsv 

	perl $codedir/overlap_intra_TAD_ChIP_Seq_peaks.pl \
        	$outdir/final_results_unchanged-TADs.tsv \
        	$outdir/diff_bind."$object2"-vs-"$object1".q100-zscore.bed \
        	$extension \
        	$outdir/unchanged-vs-diffBind-FC_intra_TAD.tsv 

	R --no-save $outdir/inactive-vs-diffBind-FC_intra_TAD.tsv \
     	   	$outdir/active-vs-diffBind-FC_intra_TAD.tsv \
		$outdir/unchanged-vs-diffBind-FC_intra_TAD.tsv \
       		$outdir/"$object2"-vs-"$object1"_diffBind_intra_TAD < $codedir/plot_boxes.r
	endif
exit


# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------

# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."


