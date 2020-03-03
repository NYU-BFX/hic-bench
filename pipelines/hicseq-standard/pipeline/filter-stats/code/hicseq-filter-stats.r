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
library(ggplot2)
library(RColorBrewer)
library(grid)

# Read arguments
args <- commandArgs(TRUE)
output <- sprintf("%s", args[1])
filenames <- strsplit(sprintf("%s", args[2]),' ')[[1]]

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

# Create data frame
ce <- ddply(total,"SAMPLE",transform, PERCENT_READS=READS/sum(as.numeric(READS))*100)
#ce <- subset(ce,ce$READS!=0)
ce = mutate(ce, READ_CATEGORY=factor(READ_CATEGORY, levels=rev(my_classes), ordered=TRUE))

# Get the output name
out <- paste(output,"percent.pdf",sep="/")

pdf(sprintf("%s",out))
ggplot(ce, aes(x=SAMPLE,y=PERCENT_READS, fill=READ_CATEGORY))+
	geom_bar(stat="identity") +
	theme(legend.background = element_rect(fill="gray90", size=.5, linetype="dotted")) + theme(legend.key.size = unit(0.4,"cm"), plot.margin = unit(c(3.5,0.5,0.5,0.5),"cm")) + theme(legend.position=c(0.5, 1.15)) +
	guides(fill=guide_legend(reverse=TRUE, title=NULL, ncol=4)) +
	scale_fill_manual(values=my_colors) +
	xlab("Sample") +
	theme(axis.title.x = element_text(size = rel(1.0), angle = 00)) + ylab("Reads (percentage)") + 
	theme(axis.title.y = element_text(size = rel(1.0), angle = 90)) +
	theme(axis.text.x = element_text(size=9, angle = 60, hjust = 1)) + 
	scale_y_continuous(breaks=seq(0,100,10), expand=c(0,0)) + coord_flip()
dev.off()


# Repeat for raw reads
ce1 <- ddply(total,"SAMPLE",transform, READS)
ce1 <- subset(ce,ce$READS!=0)

# Save to file
out <- paste(output,"counts.pdf",sep="/")

pdf(sprintf("%s",out))
ggplot(ce1, aes(x=SAMPLE,y=READS/1000000, fill=READ_CATEGORY))+
	geom_bar(stat="identity") +
	theme(legend.background = element_rect(fill="gray90", size=.5, linetype="dotted")) + theme(legend.key.size = unit(0.4,"cm"), plot.margin = unit(c(3.5,0.5,0.5,0.5),"cm")) + theme(legend.position=c(0.5, 1.15)) +
	guides(fill=guide_legend(reverse=TRUE, title=NULL, ncol=4)) +
	scale_fill_manual(values=my_colors) +
	xlab("Sample") +
	theme(axis.title.x = element_text(size = rel(1.0), angle = 00)) +
        ylab("Reads (million)") + 
	theme(axis.title.y = element_text(size = rel(1.0), angle = 90)) +
	theme(axis.text.x = element_text(size=9, angle = 60, hjust = 1)) + 
	scale_y_continuous(expand=c(0,0)) + coord_flip()
dev.off()






