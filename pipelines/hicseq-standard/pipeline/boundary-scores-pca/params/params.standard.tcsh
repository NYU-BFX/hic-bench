#!/bin/tcsh

source ./inputs/params/params.tcsh

set chrom_excluded = 'chr[MYX]'              # excluded chromosomes

set group_var = 'cell-type'                  # grouping variable (from sample sheet) to be used for color assignment)

set pca_params = '--show-text --use-short-names --plain'

#set label_fields = '1-'                     # indicate which dash-separated fields of the sample label will be printed on the PCA (e.g. '1-3', '1,2,4', '2-')             
set label_fields = '1-2'                     # indicate which dash-separated fields of the sample label will be printed on the PCA (e.g. '1-3', '1,2,4', '2-')             

