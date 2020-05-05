#!/usr/bin/Rscript
#$ -S /usr/bin/Rscript

suppressPackageStartupMessages(library(optparse))

usage = "\
\tRscript virtual4C.r [OPTIONS] OUTPUT-DIRECTORY CHROMOSOME SPARSE-MATRIX"

##
## Rscript ./code/virtual4C.r --vp-file='protein_coding2.bed' out chr8 X/matrix.chr8.mtx
##

option_list <- list(
  make_option(c("--vp-file"),default="protein_coding.bed", help="viewpoint bed file"),
  make_option(c("--nreads"),default=0, help="number of sequenced read pairs in input sample"),
  make_option(c("-u","--unit"),default=0, help="maximum resolution (bp)"),
  make_option(c("-d","--maxdist"),default=2500000, help="maximum distance from viewpoint (bp)"),
  make_option(c("-r","--radius"),default=10000, help="radius around viewpoints and target anchors (bp)")
)

# process command line arguments
arguments = parse_args(args=commandArgs(trailingOnly=T), OptionParser(usage=usage,option_list=option_list), positional_arguments=c(0,Inf))
opt = arguments$options
inputs = arguments$args
if (length(inputs) != 3) { write("Error: wrong number of inputs! Use --help to see help information", stderr()); quit(save='no') }

# input parameters
outdir = inputs[1]
chrname = inputs[2]
n_reads = as.integer(opt$"nreads")      # number of sequenced read pairs in sample 1
U = as.integer(opt$"unit")              # maximum resolution (bp)
d = as.integer(opt$"maxdist")           # maximum distance from viewpoint
r = as.integer(opt$"radius")            # radius around viewpoint

# input matrices
mat1 = inputs[3]         # e.g. DP/matrix.chr8.mtx

# load libraries
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(zoo))

# adjust by unit 
R = r %/% U
D = d %/% U
W = 2*R + 1

# divide by this constant below to obtain CPK2B values (counts per kilobase^2 per billion read pairs)
CPK2B = (W*U/1000)^2*(n_reads/1000000000)

# load matrix
write("Loading matrix...",stderr())
X <- readMM(mat1)
write(format(object.size(X),units="auto",standard="SI"),file=stderr())

# make symmetric
X = X+t(X); diag(X) = diag(X)/2
Xn = nrow(X)

# set column labels by coordinate start)
colnames(X) = seq(0,ncol(X)-1)*U

# load viewpoint information
write("Loading viewpoint information...",stderr())
G = read.table(opt$"vp-file")
colnames(G) = c("chr","start","end","label","score","strand")
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

# process each viewpoint separately
write(paste("Testing",length(vp_list),"viewpoints..."),stderr())
options(scipen=999)                                                       # disable scientific notation
n_vp = length(vp_list)
for (k in 1:n_vp) 
{
  write(paste0(chrname," ",round(100*k/n_vp,2),"%"),stderr())

  # generate raw virtual 4Cs
  VP = vp_list[k]
  x = round(v4C(X,VP,R,D,W) / CPK2B,4)
  
  # generate coordinates
  coord_start = as.numeric(names(x))
  coord_end = coord_start + U
  
  # generate v4C files
  x_out = cbind(coord_start,coord_end,x)  
  x_out = x_out[!is.na(x_out[,3]),]
  x_out = cbind(chrname,x_out)
  filename = paste(outdir,'/',G$label[k],'-',chrname,'-v4C.bedgraph',sep='') 
  cat(paste("track type=bedGraph name=",G$label[k],"-",chrname,"\n",sep=""),file=filename)
  write.table(file=filename,x_out,append=T,quote=F,col.names=F,row.names=F,sep='\t') 
}

write("Done.",stderr())

quit(save='no')

