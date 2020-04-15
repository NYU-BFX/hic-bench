#!/usr/bin/Rscript
#$ -S /usr/bin/Rscript

suppressPackageStartupMessages(library(optparse))

usage = "\
\tRscript virtual4C.r [OPTIONS] OUTPUT-DIRECTORY CHROMOSOME SPARSE-MATRIX"

##
## Rscript ./code/virtual4C.r --gene-file='protein_coding2.bed' out chr8 X/matrix.chr8.mtx
##

option_list <- list(
  make_option(c("--gene-file"),default="protein_coding.bed", help="gene bed file"),
  make_option(c("-s","--scale"),default=1.0, help="scaling factor"),
  make_option(c("-d","--maxdist"),default=2500000, help="maximum distance from viewpoint (bp)"),
  make_option(c("-r","--radius"),default=10000, help="radius around viewpoint (bp)"),
  make_option(c("-w","--window"),default=20000, help="size of rolling window (bp)")
)

# process command line arguments
arguments = parse_args(args=commandArgs(trailingOnly=T), OptionParser(usage=usage,option_list=option_list), positional_arguments=c(0,Inf))
opt = arguments$options
inputs = arguments$args
if (length(inputs) != 3) { write("Error: wrong number of inputs! Use --help to see help information", stderr()); quit(save='no') }

# input parameters
outdir = inputs[1]
chrname = inputs[2]
scaling = as.numeric(opt$"scale")       # scaling factor (for sequencing depth adjustment)
d = as.integer(opt$"maxdist")           # maximum distance from viewpoint
r = as.integer(opt$"radius")            # radius around viewpoint
w = as.integer(opt$"window")            # rolling window size

# input matrices
mat1 = inputs[3]         # e.g. DP/matrix.chr8.mtx

suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(zoo))

# adjust by unit 
U = 100            # 100bp is the finest resolution possible: this is hard-coded in the matrix-sparse step
R = r %/% U
D = d %/% U
W = w %/% U

# load matrix 1
write("Loading matrix...",stderr())
X <- readMM(mat1)
X = X+t(X); diag(X) = diag(X)/2
Xn = nrow(X)
write(format(object.size(X),units="auto",standard="SI"),file=stderr())

# load viewpoint information
write("Loading viewpoint information...",stderr())
G = read.table(opt$"gene-file")
colnames(G) = c("chr","start","end","gene","score","strand")
G = G[G$chr==chrname,,drop=FALSE]
if (nrow(G)==0) { write(paste("No viewpoints found on chromosome ",chrname,".",sep=''),stderr()); quit(save='no') }

# virtual 4C functions
v4C = function(X,VP,R,D,W) 
{
  x = colSums(X[max(VP-R,1):min(VP+R,Xn),max(VP-D,1):min(VP+D,Xn)])
  x = rollsum(x, k=W, align="center", fill=NA)
  return (x)
}

# determine viewpoints
vp_list = as.numeric(apply(G,1,function(v) { if (v["strand"]=='+') { v["start"] } else { v["end"] } })) %/% U

write(paste("Testing",length(vp_list),"viewpoints..."),stderr())
options(scipen=999)                                                       # disable scientific notation
for (k in 1:length(vp_list)) 
{
  # generate raw virtual 4Cs
  VP = vp_list[k]
  x = round(scaling * v4C(X,VP,R,D,W),4)
  
  # generate coordinates
  coord_end = U*(max(VP-D,1):min(VP+D,Xn))
  coord_start = coord_end - U
  
  # generate v4C files
  x_out = cbind(coord_start,coord_end,x)  
  x_out = x_out[!is.na(x_out[,3]),]
  x_out = cbind(chrname,x_out)
  filename = paste(outdir,'/',G$gene[k],'-',chrname,'-v4C.bedgraph',sep='') 
  cat(paste("track type=bedGraph name=",G$gene[k],"-",chrname,"\n",sep=""),file=filename)
  write.table(file=filename,x_out,append=T,quote=F,col.names=F,row.names=F,sep='\t') 
}

write("Done.",stderr())

quit(save='no')

