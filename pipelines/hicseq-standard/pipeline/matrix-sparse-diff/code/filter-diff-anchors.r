#!/usr/bin/Rscript
#$ -S /usr/bin/Rscript

suppressPackageStartupMessages(library(optparse))

usage = "\
\tRscript filter-diff-anchors.r [OPTIONS] DIFF-ANCHORS-CSV-FILE"

##
##  Rscript filter-diff-anchors.r --min-dist=10000 --min-val=4.0 diff-anchors.csv >! filtered-diff-anchors.csv
##

######################################################################
## 
## MAIN
## 
######################################################################
option_list <- list(
  make_option(c("-o","--output-dir"),default=".", help="output directory"),
  make_option(c("--min-dist"),default=10000, help="minimum anchor-pair distance (bp)"),
  make_option(c("--min-val"),default=4.0, help="minimum anchor-pair value (CPK2B: Counts Per Kilobase^2 per Billion read pairs)")
)

# process command line arguments
arguments = parse_args(args=commandArgs(trailingOnly=T), OptionParser(usage=usage,option_list=option_list), positional_arguments=c(0,Inf))
opt = arguments$options
inputs = arguments$args
if (length(inputs) != 1) { write("Error: wrong number of inputs! Use --help to see help information", stderr()); quit(save='no') }

# input parameters
outdir = opt$"output-dir"                           # output directory
min_dist = as.numeric(opt$"min-dist")               # minimum anchor-pair distance
min_val = as.integer(opt$"min-val")                 # minimum anchor-pair CPK2B value
pseudo = 1.0

# input data
write("Reading input anchor file...",stderr())
inp_file = inputs[1]
inp_data = read.csv(inp_file,header=T,check.names=FALSE)
write("Removing bias and filtering anchor pairs...",stderr())

# calculate bias
flt_dist = abs(inp_data[,"Anchor distance"])>=min_dist
mean_values = apply(inp_data[flt_dist,c("CPK2B sample 1","CPK2B sample 2")],2,mean)
bias = mean_values[2]/mean_values[1]
write(paste0("Bias = ",bias),stderr())

# filter data (uncorrected)
flt_val = (inp_data[,"CPK2B sample 1"]>=min_val) | (inp_data[,"CPK2B sample 2"]>=min_val)
filtered_uncorrected = inp_data[flt_dist&flt_val,]
filtered_uncorrected[,"Log2FC"] = log2( sapply(filtered_uncorrected[,"CPK2B sample 2"],max,pseudo) / sapply(filtered_uncorrected[,"CPK2B sample 1"],max,pseudo) )

# correct the bias and filter
filtered_corrected = inp_data 
filtered_corrected[,"CPK2B sample 1"] = filtered_corrected[,"CPK2B sample 1"]*bias
filtered_corrected[,"Count sample 1"] = round(filtered_corrected[,"Count sample 1"]*bias,0)
flt_val = (filtered_corrected[,"CPK2B sample 1"]>=min_val) | (filtered_corrected[,"CPK2B sample 2"]>=min_val)
filtered_corrected = filtered_corrected[flt_dist&flt_val,]
filtered_corrected[,"Log2FC"] = log2( sapply(filtered_corrected[,"CPK2B sample 2"],max,pseudo) / sapply(filtered_corrected[,"CPK2B sample 1"],max,pseudo) )

# print biases
write(paste0("Average log2 FC (uncorrected) = ",mean(filtered_uncorrected[,"Log2FC"])),stderr())
write(paste0("Average log2 FC (corrected) = ",mean(filtered_corrected[,"Log2FC"])),stderr())

# save outputs
store = function(X,fileX)
{
  X[,"CPK2B sample 1"] = round(X[,"CPK2B sample 1"],3)
  X[,"CPK2B sample 2"] = round(X[,"CPK2B sample 2"],3)
  X[,"Log2FC"] = round(X[,"Log2FC"],3)
  options(scipen=999)
  write.table(X,file=fileX,quote=F,col.names=TRUE,row.names=FALSE,sep=',')  
}
store(filtered_uncorrected,paste0(outdir,"/","filtered-diff-anchors.csv"))
store(filtered_corrected,paste0(outdir,"/","filtered-corrected-diff-anchors.csv"))

warnings()

write("Done.",stderr())

quit(save='no')

