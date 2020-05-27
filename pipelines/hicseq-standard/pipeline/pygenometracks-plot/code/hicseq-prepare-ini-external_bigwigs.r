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
bigwig_file_1 = args[1]
assayobject = args[2]
color_1 = args[3]
index = as.numeric(args[4])
outdir = args[5]

params.template = readLines("params/params.template.for.external_bigwigs.txt")

bigwig_file_index = grep("template_and_modify_bigwig_file", params.template)
params.template[bigwig_file_index] = gsub("template_and_modify_bigwig_file", bigwig_file_1, params.template[bigwig_file_index])

bigwig_title_index = grep("template_and_modify_bigwig_title", params.template)
params.template[bigwig_title_index] = gsub("template_and_modify_bigwig_title", assayobject, params.template[bigwig_title_index])

bigwig_color_index = grep("template_and_modify_bigwig_color", params.template)
params.template[bigwig_color_index] = gsub("template_and_modify_bigwig_color", color_1, params.template[bigwig_color_index])

writeLines(params.template, con=paste0(outdir, "/", sprintf("%02d", index), ".", assayobject, "-bigwig.made.ini"))







