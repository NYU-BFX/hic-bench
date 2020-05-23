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
library(plyr)
library(grid)

# Read arguments
args <- commandArgs(TRUE)
filenames = dir(path=args[1], full.names=T)
object = args[2]

# Now create a datalist where all the stats
# will be stored

data <- list()

for (i in 1:length(filenames)) {
	df <- data.frame(read.table(sprintf("%s/stats.tsv", filenames[i]), header=TRUE, stringsAsFactors=FALSE))
    #df$i <- i  # maybe you want to keep track of which iteration produced it?
	df$SAMPLE <- rep(as.character(basename(filenames[i])),dim(df)[1])
	colnames(df) <- c("READ_CATEGORY","READS","PERCENT","SAMPLE")
	data[[i]] <- df
}

# Get all the data
total <- do.call(rbind, data)
total$SAMPLE <- factor(total$SAMPLE)

# Find number of unique categories
read_category_no <- length(unique(total$READ_CATEGORY))

# Assign color palette based on category number
my_classes = c("ds-accepted-intra","ds-accepted-inter","ds-duplicate-intra","ds-duplicate-inter","multihit","single-sided","ds-no-fragment","ds-same-fragment","ds-too-close","ds-too-far","unpaired","unmapped","unclassified")
my_colors <- c("#33a02c","#b2df8a","#e31a1c","#fb9a99","#a6cee3","#1f78b4","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99","#b15928","#000000")
names(my_colors) = my_classes

# Create data frame for raw reads
ce1 <- ddply(total,"SAMPLE",transform, READS)
# Remove "rep1 or rep2" texts
ce1$group = gsub("rep.*$", "", ce1$SAMPLE)
# Remove Arima, DpnII, HindIII, MboI, NcoI texts
ce1$group = gsub("(-Arima-)|(-DpnII-)|(-HindIII-)|(-MboI-)|(-NcoI-)", "", ce1$group)

# Only need ds-accepted-intra.
output.needed = subset(ce1, READ_CATEGORY == "ds-accepted-intra" & group == object)$READS

# Print the scaling factor to produce results in 1000 * CPM.
options(digits=15)
scale.factor = 1 / sum(output.needed) * 1e9

cat(scale.factor)





