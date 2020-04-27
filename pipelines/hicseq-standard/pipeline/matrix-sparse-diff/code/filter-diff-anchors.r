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
  make_option(c("--min-dist"),default=10000, help="minimum anchor-pair distance (bp)"),
  make_option(c("--min-val"),default=4.0, help="minimum anchor-pair value (CPK2B: Counts Per Kilobase^2 per Billion read pairs)")
)

# process command line arguments
arguments = parse_args(args=commandArgs(trailingOnly=T), OptionParser(usage=usage,option_list=option_list), positional_arguments=c(0,Inf))
opt = arguments$options
inputs = arguments$args
if (length(inputs) != 1) { write("Error: wrong number of inputs! Use --help to see help information", stderr()); quit(save='no') }

# input parameters
min_dist = as.numeric(opt$"min-dist")               # minimum anchor-pair distance
min_val = as.integer(opt$"min-val")                 # minimum anchor-pair CPK2B value
correct_bias = FALSE

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

# calculate mean log2 fold-changes on uncorrected filtered data
flt_val = (abs(inp_data[,"CPK2B sample 1"])>=min_val) | (abs(inp_data[,"CPK2B sample 2"])>=min_val)
m = mean(inp_data[flt_dist&flt_val,"Log2FC"])
write(paste0("Average log2 FC (uncorrected) = ",m),stderr())

if (correct_bias == FALSE) {
  filtered_data = inp_data[flt_dist&flt_val,]

} else {
  # correct the bias and filter
  inp_data[,"CPK2B sample 1"] = inp_data[,"CPK2B sample 1"]*bias
  flt_val = (abs(inp_data[,"CPK2B sample 1"])>=min_val) | (abs(inp_data[,"CPK2B sample 2"])>=min_val)
  filtered_data = inp_data[flt_dist&flt_val,]

  # correct the log2 fold-changes
  pseudo = 0.5
  filtered_data[,"Log2FC"] = log2( sapply(filtered_data[,"CPK2B sample 2"],max,pseudo)/sapply(filtered_data[,"CPK2B sample 1"],max,pseudo))
  m = mean(filtered_data[,"Log2FC"])
  write(paste0("Average log2 FC (corrected) = ",m),stderr())
}

# format table
filtered_data[,"CPK2B sample 1"] = round(filtered_data[,"CPK2B sample 1"],3)
filtered_data[,"CPK2B sample 2"] = round(filtered_data[,"CPK2B sample 2"],3)
filtered_data[,"Log2FC"] = round(filtered_data[,"Log2FC"],3)

# write output
options(scipen=999)
write.table(filtered_data,quote=F,col.names=TRUE,row.names=FALSE,sep=',')  

warnings()

write("Done.",stderr())

quit(save='no')

