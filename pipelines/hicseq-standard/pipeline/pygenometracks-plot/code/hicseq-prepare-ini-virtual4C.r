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
virtual4C_file_1 = args[1]
virtual4C_color_1 = args[2]
index = as.numeric(args[3])
object = args[4]
virtual4C_overlay = args[5]
outdir = args[6]

params.template = readLines("params/params.template.for.virtual4C.txt")

virtual4C_file_index = grep("template_and_modify_virtual4C_file", params.template)
params.template[virtual4C_file_index] = gsub("template_and_modify_virtual4C_file", virtual4C_file_1, params.template[virtual4C_file_index])

virtual4C_color_index = grep("template_and_modify_virtual4C_color", params.template)
params.template[virtual4C_color_index] = gsub("template_and_modify_virtual4C_color", virtual4C_color_1, params.template[virtual4C_color_index])

virtual4C_overlay_index = grep("template_and_modify_virtual4C_overlay", params.template)
params.template[virtual4C_overlay_index] = gsub("template_and_modify_virtual4C_overlay", virtual4C_overlay, params.template[virtual4C_overlay_index])


writeLines(params.template, con=paste0(outdir, "/", sprintf("%02d", index), ".virtual4C.", object, ".made.ini"))







