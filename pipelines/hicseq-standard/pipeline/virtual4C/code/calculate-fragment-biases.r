#!/usr/bin/Rscript
#$ -S /usr/bin/Rscript

suppressPackageStartupMessages(library(optparse))

usage = "\
\tcalculate-fragment-biases.r [OPTIONS] OUTPUT-FILE SPARSE-MATRIX"

##
## Rscript calculate-fragment-biases.r out X/matrix.chr8.mtx
##

option_list <- list(
  make_option(c("-w","--window"),default=20000, help="size of rolling window (bp)")
)

# process command line arguments
arguments = parse_args(args=commandArgs(trailingOnly=T), OptionParser(usage=usage,option_list=option_list), positional_arguments=c(0,Inf))
opt = arguments$options
inputs = arguments$args
if (length(inputs) != 2) { write("Error: wrong number of inputs! Use --help to see help information", stderr()); quit(save='no') }

# input parameters
outfile = inputs[1]
mat1 = inputs[2]         # e.g. DP/matrix.chr8.mtx
w = as.integer(opt$"window")            # rolling window size

# load libraries
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(zoo))

# adjust by unit 
U = 100            # 100bp is the finest resolution possible: this is hard-coded in the matrix-sparse step
W = w %/% U

# load matrix
write("Loading matrix...",stderr())
x <- readMM(mat1)
n = ncol(x)

# determine chromosome name
chrname = 'chr10'

# calculate rolling sum
write("Computing rolling sum...",stderr())
x = rollsum(x[1,], k=W, align="center", fill=NA)
  
# generate coordinates
options(scipen=999)                                                       # disable scientific notation
coord_end = U*(1:n)
coord_start = coord_end - U
  
# generate v4C files
x_out = cbind(coord_start,coord_end,x)  
x_out = x_out[!is.na(x_out[,3]),]
x_out = cbind(chrname,x_out)

# store as bedgraph file
cat("track type=bedGraph\n",file=outfile)
write.table(file=outfile,x_out,append=T,quote=F,col.names=F,row.names=F,sep='\t') 

write("Done.",stderr())

quit(save='no')

