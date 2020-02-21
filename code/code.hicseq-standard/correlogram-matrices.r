#!/bin/Rscript

# This is a script that takes as input the table with the sample pairs, the distances, lambdas
# and correlations and outputs the correlograms for each distance and each lambda, as well as
# the summary tables for each distance

# Check for required packages
# and install
for (package in c("plyr","reshape2","RColorBrewer","corrplot")) {
        if(package %in% rownames(installed.packages()) == FALSE){install.packages(package, repos="http://cran.us.r-project.org")}
}

# Load the required packages
library(plyr)
library(reshape2)
library(corrplot)

args <- commandArgs(trailingOnly = TRUE)

if(length(args)!=1){print("USAGE: Rscript correlogram-matrices.r INPUT-TABLE"); quit(save="no")}


# read the data
data <- read.table(sprintf("%s", args[1]), header=TRUE, sep="\t", check.names=FALSE)

# Find the number of columns
column_no <- length(colnames(data))

# Get the part that is standard
data1 <- data[,1:6]

# Get the part with d
data2 <- data.frame(data[,7:column_no])

# Get the number of d
d_no <- length(colnames(data2))
d_names <- colnames(data)[7:column_no]
colnames(data2) <- d_names

# Get the name of directory
dir <- dirname(args[1])

# Now combine the standard with the variable part
out2 <- paste(dir,"correlograms.pdf",sep="/")
pdf(out2)
for (i in 1:d_no) {
	data3 <- data.frame(cbind(data1,data2[,i]))
	colnames(data3) <- c("sample1","sample2","pair","method","chrom","lambda","correlation")
	d1 <- ddply(data3, c("sample1", "sample2", "pair","method","lambda"), summarise, mean=mean(correlation), min=min(correlation), max=max(correlation), 
	N = length(correlation), sd = sd(correlation), se = sd / sqrt(N))
	#Export the file for each one of the distances
	out <- paste(d_names[i],"summary.tsv",sep="-")
	out1 <- paste(dir,out,sep="/")
	write.table(d1, file=out1, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
	#Get the unique lambdas
	lambdas <- unique(d1$lambda)
	for (j in 1:length(lambdas)) {
        	d2 <- subset(d1, d1$lambda==lambdas[j])
  		m <- acast(d2, sample1~sample2, value.var='mean')
  		col <- colorRampPalette(c("#4477AA", "#77AADD", "#FFFFFF", "#EE9988", "#BB4444"))
		title_part1 <- paste("distance", d_names[i], sep="=")
        	title_part2 <- paste("lambda",lambdas[j],sep="=")
		title <- paste(title_part1, title_part2, sep=",")
		cat("Creating correlogram for", title_part1, title_part2, "...\n")
  		corrplot(m, method="color", color=col(200), type="full", title=sprintf("%s", title), order="original", addCoef.col = "black", tl.cex=0.75, tl.col="black", tl.srt=45, mar=c(0,0,1,0), diag=FALSE)
	}
}
dev.off()
