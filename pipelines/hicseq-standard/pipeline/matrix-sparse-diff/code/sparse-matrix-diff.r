#!/usr/bin/Rscript
#$ -S /usr/bin/Rscript

suppressPackageStartupMessages(library(optparse))

usage = "\
\tRscript sparse-matrix-diff.r [OPTIONS] OUTPUT_DIRECTORY CHROMOSOME SPARSE-MATRIX-1 SPARSE-MATRIX-2 LABEL-1 LABEL-2"

##
## Rscript ./code/sparse-matrix-diff.r --gene-file='protein_coding2.bed' out chr8 X/matrix.chr8.mtx Y/matrix.chr8.mtx DP TALL
##

######################################################################
## 
## FUNCTION findSegments | Find segments above a threshold
## 
######################################################################
findSegments <- function(x, y, min_diff, min_region_size=5000%/%U) 
{
  
  findSegmentsAux <- function(dxy, threshold) 
  {
    hit <- which(dxy > threshold)
    n <- length(hit)
    ind <- which(hit[-1] - hit[-n] > 1)
	if (length(ind)==0) return(c())
    starts <- c(hit[1], hit[ ind+1 ])
    ends <- c(hit[ ind ], hit[n])
    size <- ends - starts + 1
    return(cbind(starts,ends,size))
  }

  # init
  dxy = y - x
  seg = seg1 = seg2 = c()
  
  # find maxima
  dd = rollmax(dxy,k=min_region_size,na.pad=T,align="center")
  if (max(dd,na.rm=T)>=min_diff) seg1 = findSegmentsAux(dd,min_diff)

  # find minima
  dd = rollmax(-dxy,k=min_region_size,na.pad=T,align="center")
  if (max(dd,na.rm=T)>=min_diff) seg2 = rbind(seg, findSegmentsAux(dd,min_diff))

  #print("diagnostics:");  print(seg1);   print(seg2)
  
  seg = rbind(seg1,seg2)
  if (is.null(seg)==FALSE) {
    xvalue = apply(seg,1,function(v) mean(x[v[1]:v[2]]))           # add mean value for sample x
    yvalue = apply(seg,1,function(v) mean(y[v[1]:v[2]]))           # add mean value for sample y
	dvalue = yvalue - xvalue                                       # add difference of mean values
	seg = cbind(seg,xvalue,yvalue,dvalue)
  }
 
  return(seg)
}


######################################################################
## 
## MAIN
## 
######################################################################
option_list <- list(
  make_option(c("--gene-file"),default="protein_coding.bed", help="gene bed file"),
  make_option(c("-d","--maxdist"),default=2500000, help="maximum distance from viewpoint (bp)"),
  make_option(c("-w","--window"),default=20000, help="size of rolling window (bp)"),
  make_option(c("--mincount"),default=50, help="minimum viewpoint count for virtual 4Cs"),
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
d = as.integer(opt$maxdist)           # maximum distance from viewpoint
w = as.integer(opt$window)            # rolling window size
mincount = as.integer(opt$mincount)   # minimum viewpoint count for virtual 4Cs
mindiff = as.numeric(opt$mindiff)     # minimum difference
r = 5000                              # radius around viewpoint (bp)

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
Xn = nrow(X)
write(format(object.size(X),units="auto",standard="SI"),file=stderr())

# load matrix 2
write("Loading matrix 2...",stderr())
Y <- readMM(mat2)
Y = Y+t(Y); diag(Y) = diag(Y)/2
Yn = nrow(Y)
write(format(object.size(Y),units="auto",standard="SI"),file=stderr())

# check matrix sizes
if (Xn!=Yn) { write("Error: input matrices have different sizes!\n",file=stderr()); quit(save='no') }

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

# identify differential regions
regdiff = function(z,min_diff,min_region_size=1000%/%U)
{
  z[z<min_diff] = NA
  sum(diff(which(!is.na(z)))>=min_region_size)
}


# determine viewpoints
vp_list = as.numeric(apply(G,1,function(v) { if (v["strand"]=='+') { v["start"] } else { v["end"] } })) %/% U

# initialize viewpoint stats matrix
col_labels = c( 
	paste(L1,"max count"), 
	paste(L2,"max adjusted count"), 
	paste(L1,"max2 count"), 
	paste(L2,"max2 adjusted count"), 
	paste(L1,"sum counts"), 
	paste(L2,"sum adjusted counts"), 
	paste(L1,"-high regions (scaled1)",sep=''), 
	paste(L2,"-high regions (scaled1)",sep=''), 
	paste(L1,"-high regions (scaled2)",sep=''), 
	paste(L2,"-high regions (scaled2)",sep=''), 
	"Max count", 
	"Max2 count", 
	"Log2FC",
	"Diff area (scaled1)",
	"Diff area (scaled2)",
	"Absolute diff area (scaled1)",
	"Absolute diff area (scaled2)"
)
vp_stats = matrix(0,length(vp_list),length(col_labels))
rownames(vp_stats) = G$gene
colnames(vp_stats) = col_labels

write(paste("Testing",length(vp_list),"viewpoints..."),stderr())
file_diff_regions = paste(outdir,'/diff-regions.csv',sep='')
cat("vp-name,vp-chr,vp-start,vp-end,target-distance,target-start,target-end,target-length,sample1-value,sample2-value,diff\n",file=file_diff_regions)
for (k in 1:length(vp_list)) 
{
  # generate raw virtual 4Cs
  VP = vp_list[k]
  vp_label = as.character(G$gene[k])
  x = v4C(X,VP,R,D,W)
  y = v4C(Y,VP,R,D,W)
  
  # generate coordinates
  coord_end = U*(max(VP-D,1):min(VP+D,Xn))
  coord_start = coord_end - U
  
  # generate max & sum stats
  x_max = max(x)                                                    # max value in sample X
  y_max = max(y)                                                    # max value in Y (adjusted)
  imax = which.max(x)                                               # max value outside the main peak area in sample X
  xx = x; xx[max(1,imax-20*R):min(imax+20*R,Xn)] = NA
  x_max2 = max(xx,na.rm=T)
  imax = which.max(y)                                               # max value outside the main peak area in sample Y
  yy = y; yy[max(1,imax-20*R):min(imax+20*R,Xn)] = NA
  y_max2 = max(yy,na.rm=T)
  x_sum = sum(x)                                                    # sum values in sample X
  y_sum = sum(y)                                                    # sum values in Y (adjusted)
  max_xy = max(x_max,y_max)                                         # max value across X and Y
  max2_xy = max(x_max2,y_max2)                                      # max2 value across X and Y

  # scale to max and maxmax and calculate differences
  xs1 = x/x_max
  ys1 = y/y_max
  dxy1 = ys1 - xs1
  xs2 = x/max_xy
  ys2 = y/max_xy
  dxy2 = ys2 - xs2
  
  # generate diff regions stats
  nreg_dx1 = regdiff(-dxy1,mindiff)                                 # max-scaling: number of regions higher in sample X (above tolerance)
  nreg_dy1 = regdiff(dxy1,mindiff)                                  # max-scaling: number of regions higher Y
  nreg_dx2 = regdiff(-dxy2,mindiff)                                 # maxmax-scaling: number of regions higher in sample X (above tolerance)
  nreg_dy2 = regdiff(dxy2,mindiff)                                  # maxmax-scaling: number of regions higher Y

  # generate diff metrics
  logfc = log2(y_sum/x_sum)                                         # log2 ratio (Y vs X)
  delta_area_scaled1 = sum(dxy1)/max(sum(ys1),sum(xs1))             # normalized total max-scaled difference (Y vs X)
  delta_area_scaled2 = sum(dxy2)/max(sum(ys2),sum(xs2))             # normalized total maxmax-scaled difference (Y vs X)
  abs_delta_area_scaled1 = sum(abs(dxy1))/max(sum(ys1),sum(xs1))    # normalized total absolute max-scaled difference (Y vs X)
  abs_delta_area_scaled2 = sum(abs(dxy2))/max(sum(ys2),sum(xs2))    # normalized total absolute maxmax-scaled difference (Y vs X)

  # generate differential regions
  seg = findSegments(xs2,ys2,min_diff=mindiff)                      # use max-max scaling for differential regions
  if (is.null(seg)==FALSE) {
    vp_chr = chrname
	vp_start = (VP-1)*U
	vp_end = (VP+1)*U
    seg[,1] = coord_start[seg[,1]]                                  # fix coordinates
    seg[,2] = coord_end[seg[,2]]
	seg[,3] = seg[,2]-seg[,1]
    loop_dist = (seg[,1]+seg[,2])/2-mean(vp_start,vp_end)           # target anchor distance from viewpoint
	seg = cbind(vp_label,vp_chr,vp_start,vp_end,loop_dist,seg)
    write.table(file=file_diff_regions,seg,append=T,quote=F,col.names=F,row.names=F,sep=',') 
  }
 
  # generate v4C files
  if (max_xy>=mincount) {
    filename = paste(outdir,'/',vp_label,'-',chrname,'-v4C.csv',sep='') 
    dataset = round(cbind(x,y,xs1,ys1,dxy1),3)
    colnames(dataset) = c( paste(L1,"counts"), paste(L2,"adjusted counts"), paste(L1,"max-scaled"), paste(L2,"max-scaled"), "Diff-scaled") 
    write.table(file=filename,dataset,quote=F,col.names=T,row.names=F,sep=',') 
  }
  
  # store stats 
  vp_stats[k,] = c( 
		x_max, 
		round(y_max,0), 
		x_max2, 
		round(y_max2,0), 
		x_sum, 
		round(y_sum,0), 
		nreg_dx1, 
		nreg_dy1, 
		nreg_dx2, 
		nreg_dy2, 
		round(max_xy,0), 
		round(max2_xy,0), 
		round(logfc,3), 
		round(delta_area_scaled1,3), 
		round(delta_area_scaled2,3), 
		round(abs_delta_area_scaled1,3), 
		round(abs_delta_area_scaled2,3)
	)
}

# write output
write.table(file=paste(outdir,'/stats.csv',sep=''),vp_stats,quote=F,col.names=NA,row.names=T,sep=',')                                 # NOTE: add GeneName to column labels

write("Done.",stderr())

quit(save='no')

