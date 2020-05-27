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
domains_diff_file_1 = args[1]
domains_diff_title_1 = args[2]
outdir = args[3]

params.template = readLines("params/params.template.for.domains_diff.txt")

domains_diff_file_index = grep("template_and_modify_domains_diff_file", params.template)
params.template[domains_diff_file_index] = gsub("template_and_modify_domains_diff_file", domains_diff_file_1, params.template[domains_diff_file_index])

domains_diff_title_index = grep("template_and_modify_domains_diff_title", params.template)
params.template[domains_diff_title_index] = gsub("template_and_modify_domains_diff_title", domains_diff_title_1, params.template[domains_diff_title_index])

writeLines(params.template, con=paste0(outdir, "/", domains_diff_title_1, ".domains_diff.made.ini"))







