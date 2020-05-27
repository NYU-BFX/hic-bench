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





