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
compartment_file_1 = args[1]
compartment_title_1 = args[2]
index = as.numeric(args[3])
outdir = args[4]

print(args)

params.template = readLines("params/params.template.for.compartments.txt")

compartment_file_index = grep("template_and_modify_compartments_file", params.template)
params.template[compartment_file_index] = gsub("template_and_modify_compartments_file", compartment_file_1, params.template[compartment_file_index])

compartment_title_index = grep("template_and_modify_compartments_title", params.template)
params.template[compartment_title_index] = gsub("template_and_modify_compartments_file", compartment_title_1, params.template[compartment_title_index])

writeLines(params.template, con=paste0(outdir, "/", sprintf("%02d", index), ".", compartment_title_1, "-hic_matrix.made.ini"))







