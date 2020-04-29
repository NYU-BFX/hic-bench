#!/bin/tcsh

source ./inputs/params/params.tcsh

set chrom_excluded = 'chr[MYX]'                    # excluded chromosomes
set group_var = 'group'                            # grouping variable (from sample sheet) to be used for color assignment)
set pca_params = '--show-text --use-short-names'
set label_fields = '1,2'                           # indicate which dash-separated fields of the sample label will be printed on the PCA (e.g. '1-3', '1,2,4', '2-')             
set centrotelo_file = $genome_dir/centrotelo.bed
set delta_cut = 1.2			           # relative delta cutoff for compartments switch classification

