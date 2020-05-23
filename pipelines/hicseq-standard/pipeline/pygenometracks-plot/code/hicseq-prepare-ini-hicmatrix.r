#!/bin/Rscript

########Description######################
#This is an R script that accepts
#a .csv file as an argument which
#contains reads for different samples
#and categories (mapped, non-mapped etc.)
#and generates the corresponding stacked
#boxplot in .pdf format
#########################################

# Check for required packages
# and install
#for (package in c("plyr","ggplot2","RColorBrewer","grid")) {
#	if(package %in% rownames(installed.packages()) == FALSE){install.packages(package, repos="http://cran.us.r-project.org")}
#}
 
# Load the libraries
# Read arguments
args = commandArgs(T)
hic_matrix_file_1 = args[1]
object = args[2]
scaling_factor_for_object = args[3]
domains_bed = args[4]
outdir = args[5]


params.template = readLines("params/params.template.for.hic_matrix.txt")

hic_matrix_file_index = grep("template_and_modify_hic_matrix_file", params.template)
params.template[hic_matrix_file_index] = gsub("template_and_modify_hic_matrix_file", hic_matrix_file_1, params.template[hic_matrix_file_index])

hic_matrix_title_index = grep("template_and_modify_hic_matrix_title", params.template)
params.template[hic_matrix_title_index] = gsub("template_and_modify_hic_matrix_title", object, params.template[hic_matrix_title_index])

hic_matrix_scale_factor_index = grep("template_and_modify_hic_matrix_scale_factor", params.template)
params.template[hic_matrix_scale_factor_index] = gsub("template_and_modify_hic_matrix_scale_factor", scaling_factor_for_object, params.template[hic_matrix_scale_factor_index])

domains_file_index = grep("template_and_modify_domains_file", params.template)
params.template[domains_file_index] = gsub("template_and_modify_domains_file", domains_bed, params.template[domains_file_index])

writeLines(params.template, con=paste0(outdir, "/", "hic_matrix.made.ini"))







