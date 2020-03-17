#!/usr/bin/Rscript
#$ -S /usr/bin/Rscript

suppressPackageStartupMessages(library(optparse))

usage = "\
\tRscript sparse-matrix-diff.r [OPTIONS] OUTPUT_DIRECTORY CHROMOSOME SPARSE-MATRIX-1 SPARSE-MATRIX-2 LABEL-1 LABEL-2"

##
## Rscript ./code/sparse-matrix-diff.r --gene-file='protein_coding2.bed' out chr8 X/matrix.chr8.mtx Y/matrix.chr8.mtx DP TALL
##

option_list <- list(
  make_option(c("-d","--maxdist"),default=2500000, help="maximum distance from viewpoint (bp)"),
  make_option(c("-w","--window"),default=20000, help="size of rolling window (bp)"),
  make_option(c("--gene-file"),default="protein_coding.bed", help="gene bed file"),
  make_option(c("--mindiff"),default=0.1, help="minimum difference (fraction)")
)

# process command line arguments
arguments = parse_args(args=commandArgs(trailingOnly=T), OptionParser(usage=usage,option_list=option_list), positional_arguments=c(0,Inf))
opt = arguments$options
inputs = arguments$args
if (length(inputs) != 6) { write("Error: wrong number of inputs! Use --help to see help information", stderr()); quit(save='no') }

# input parameters
outdir = inputs[1]
chrname = inputs[2]
d = as.integer(opt$maxdist)        # maximum distance from viewpoint
w = as.integer(opt$window)         # rolling window size
mindiff = as.numeric(opt$mindiff)  # minimum difference
r = 5000                           # radius around viewpoint (bp)

# input matrices
mat1 = inputs[3]         # e.g. DP/matrix.chr8.mtx
mat2 = inputs[4]         # e.g. TALL/matrix.chr8.mtx
L1 = inputs[5]           # e.g. DP
L2 = inputs[6]           # e.g. TALL

suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(zoo))

# adjust by unit 
U = 100            # 100bp is the finest resolution possible: this is hard-coded in the matrix-sparse step
R = r %/% U
D = d %/% U
W = w %/% U

# load matrix 1
write("Loading matrix 1...",stderr())
X <- readMM(mat1)
X = X+t(X); diag(X) = diag(X)/2
write(format(object.size(X),units="auto",standard="SI"),file=stderr())

# load matrix 2
write("Loading matrix 2...",stderr())
Y <- readMM(mat2)
Y = Y+t(Y); diag(Y) = diag(Y)/2
write(format(object.size(Y),units="auto",standard="SI"),file=stderr())

# adjust counts in second sample
a = sum(X)/sum(Y)
Y = a*Y

# load gene information
G = read.table(opt$"gene-file")
colnames(G) = c("chr","start","end","gene","score","strand")
G = G[G$chr==chrname,,drop=FALSE]
if (nrow(G)==0) { write(paste("No viewpoints found on chromosome ",chrname,".",sep=''),stderr()); quit(save='no') }

# virtual 4C functions
v4C = function(X,VP,R,D,W) 
{
  Xn = nrow(X)
  x = colSums(X[max(VP-R,1):min(VP+R,Xn),max(VP-D,1):min(VP+D,Xn)])
  x = rollsum(x, k=W, align="center")
  return (x)
}

# determine viewpoints
vp_list = as.numeric(apply(G,1,function(v) { if (v["strand"]=='+') { v["start"] } else { v["end"] } })) %/% U

# initialize viewpoint stats matrix
col_labels = c( 
	paste(L1,"max count"), 
	paste(L2,"max adjusted count"), 
	paste(L1,"sum counts"), 
	paste(L2,"sum adjusted counts"), 
	paste(L1,"-high bins",sep=''), 
	paste(L2,"-high bins",sep=''), 
	"Max count", 
	"Sum count", 
	"Log2FC",
	"Diff area",
	"Diff area scaled",
	"Absolute diff area",
	"Absolute diff area scaled"
)
vp_stats = matrix(0,length(vp_list),length(col_labels))
rownames(vp_stats) = G$gene
colnames(vp_stats) = col_labels

write(paste("Testing",length(vp_list),"viewpoints..."),stderr())
for (k in 1:length(vp_list)) 
{
  # generate raw virtual 4Cs
  VP = vp_list[k]
  x = v4C(X,VP,R,D,W)
  y = v4C(Y,VP,R,D,W)
  
  # scale to max
  xs = x/max(x)
  ys = y/max(y)
  
  # v4C difference
  dxy = ys - xs
  
  # generate stats
  x_max = max(x)                                                    # max value in sample X
  y_max = max(y)                                                    # max value in Y (adjuected)
  x_sum = sum(x)                                                    # sum values in sample X
  y_sum = sum(y)                                                    # sum values in Y (adjusted)
  n_dx = sum(dxy<=-mindiff)                                         # number of bins higher in sample X (above tolerance)
  n_dy = sum(dxy>=+mindiff)                                         # number of bins higher in Y
  max_xy = max(x_max,y_max)                                         # max value across X and Y
  sum_xy = max(x_sum,y_sum)                                         # max sum value across X and Y
  logfc = log2(y_sum/x_sum)                                         # log2 ratio (Y vs X)
  delta_area = (y_sum-x_sum)/max(x_sum,y_sum)                       # normalized total difference (Y vs X)
  delta_area_scaled = sum(dxy)/max(sum(ys),sum(xs))                 # normalized total scaled difference (Y vs X)
  abs_delta_area = abs(y_sum-x_sum)/max(x_sum,y_sum)                # normalized total absolute difference (Y vs X)
  abs_delta_area_scaled = sum(abs(dxy))/max(sum(ys),sum(xs))        # normalized total absolute scaled difference (Y vs X)

  # generate v4C files
  if (x_max>=50) {
    filename = paste(outdir,'/v4C-',G$gene[k],'.csv',sep='') 
    dataset = round(cbind(x,y,xs,ys,dxy),3)
    colnames(dataset) = c( paste(L1,"counts"), paste(L2,"adjusted counts"), paste(L1,"max-scaled"), paste(L2,"max-scaled"), "Diff-scaled") 
    write.table(file=filename,dataset,quote=F,col.names=T,row.names=F,sep=',') 
  }
  
  # store stats and virtual 4C data
  vp_stats[k,] = c( x_max, y_max, x_sum, y_sum, n_dx, n_dy, max_xy, sum_xy, logfc, delta_area, delta_area_scaled, abs_delta_area, abs_delta_area_scaled )
}

# write output
write.table(file=paste(outdir,'/stats.csv',sep=''),vp_stats,quote=F,col.names=NA,row.names=T,sep=',')                                 # NOTE: add GeneName to column labels

write("Done.",stderr())

quit(save='no')

