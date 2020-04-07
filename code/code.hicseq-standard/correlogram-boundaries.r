#!/bin/Rscript

# Check for required packages
# and install
for (package in c("plyr","reshape2","RColorBrewer","corrplot")) {
        if(package %in% rownames(installed.packages()) == FALSE){install.packages(package, repos="http://cran.us.r-project.org")}
}


# Load the required libraries
library(reshape2)
library(corrplot)
library(RColorBrewer)

## different color series
col1 <- colorRampPalette(c("#7F0000","red","#FF7F00","yellow","white",
"cyan", "#007FFF", "blue","#00007F"))
col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
"#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))
col3 <- colorRampPalette(c("blue", "yellow", "red"))
col4 <- colorRampPalette(c("#7F0000","red","#FF7F00","yellow","#7FFF7F",
"cyan", "#007FFF", "blue","#00007F"))
wb <- c("white","black")

# Read the arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1){print("USAGE: Rscript correlogram-boundaries.r INPUT-FILE"); quit(save="no")}

# Read input
data <- read.table(sprintf("%s", args[1]), sep="\t", header=TRUE)
 
# Percentages
m <- acast(data, SAMPLE1~SAMPLE2~KAPPA, value.var="OVERLAP")
# Actual numbers
m1 <- acast(data, SAMPLE1~SAMPLE2~KAPPA, value.var="N.COMMON")

# Find the different kappas
kappas <- unique(data$KAPPA)

# Get the outputs
dir <- dirname(args[1])

# Create output name
raw_counts <- paste(dir,"raw_comparisons.pdf",sep="/")
corrs      <- paste(dir,"correlograms.pdf",sep="/")

# Get the raw numbers
pdf(sprintf("%s",corrs))
for (kappa in kappas) { 
	plot_title = ""   # paste("Overlap (%)", kappa, sep="-")
        m2 <- m[, , kappa] 
	m2[is.na(m2)] = 0
	corrplot.mixed(m2, lower="number", upper="pie",  tl.pos = c("lt"), tl.col="black", 
	tl.cex=0.8, tl.srt=45, diag = c("l"), cl.lim=c(0,1), color=col3(1000),
	title=plot_title, mar=c(0,0,1,0))
}
dev.off()

# Get the correlograms
pdf(sprintf("%s",raw_counts))
for (kappa in kappas) {
	plot_title = paste("Number of boundaries", kappa, sep="-")
        m3 <- m1[, , kappa] 
	corrplot(m3, mar=c(0,0,1,0), type=c("full"), is.corr=FALSE,  method=c("number"), 
		color="black", tl.pos = c("lt"), tl.col="black", tl.cex=0.8, diag = TRUE, 
		cl.pos=c("n"), title=plot_title)
}
dev.off()






