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
  make_option(c("--target-file"),default="", help="target anchors bed file (optional)"),
  make_option(c("--nreads"),default=0, help="number of sequenced read pairs in input sample"),
  make_option(c("-u","--unit"),default=0, help="maximum resolution (bp)"),
  make_option(c("-d","--maxdist"),default=2500000, help="maximum distance from viewpoint (bp)"),
  make_option(c("-r","--radius"),default=10000, help="radius around viewpoints and target anchors (bp)"),
  make_option(c("--minvalue"),default=2.0, help="minimum CPK2B (counts per kilobase^2 per billion reads) applied to virtual 5C results")
)

# process command line arguments
arguments = parse_args(args=commandArgs(trailingOnly=T), OptionParser(usage=usage,option_list=option_list), positional_arguments=c(0,Inf))
opt = arguments$options
inputs = arguments$args
if (length(inputs) != 6) { write("Error: wrong number of inputs! Use --help to see help information", stderr()); quit(save='no') }

# input parameters
outdir = inputs[1]
chrname = inputs[2]
nullRmvAnchors = inputs[4]                # TRUE: remove anchors from the null distribution
v4c_bdg = inputs[5]                       # TRUE/FALSE: produce v4c bedgraphs
normalize_bdg = inputs[6]                 # normalize bedgraph using CPK2B
n_reads = as.integer(opt$"nreads")        # number of sequenced read pairs in sample 1
U = as.integer(opt$"unit")                # maximum resolution (bp)
d = as.integer(opt$"maxdist")             # maximum distance from viewpoint
r = as.integer(opt$"radius")              # radius around viewpoint
minvalue = as.numeric(opt$"minvalue")     # minimum CPK2B (counts per kilobase^2 per billion reads) applied to virtual 5C results		
min.d.DB = 20000			  # minimum distance computed in the null distribution
max.d.DB = 1000000			  # maximum distance computed in the null distribution
bin.DB = 5000				  # bin size used to compute the null distribution

# input matrices
mat1 = inputs[3]         # e.g. DP/matrix.chr8.mtx

# load libraries
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(zoo))

# adjust by unit 
R = r %/% U
D = d %/% U
W = 2*R + 1

# adjust maximum distance (for consistency)
D = D + R*2

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

# Load viewpoint information (source anchors)
write("Loading viewpoint information...",stderr())
vp_table = read.table(opt$"vp-file")
colnames(vp_table) = c("chr","start","end","label","score","strand")
vp_table = vp_table[vp_table$chr==chrname,,drop=FALSE]
if (nrow(vp_table)==0) { write(paste("No viewpoints found on chromosome ",chrname,".",sep=''),stderr()); quit(save='no') }

# Load target anchor information
if (opt$"target-file"!="") {
  anchor_table = read.table(opt$"target-file")
  colnames(anchor_table) = c("chr","start","end","label","score","strand")
  anchor_table = anchor_table[anchor_table$chr==chrname,,drop=FALSE]
  if (nrow(anchor_table)==0) { write(paste("No target anchors found on chromosome ",chrname,".",sep=''),stderr()); quit(save='no') }
} else {
  anchor_table = vp_table
}

# Virtual 4C functions
v4C = function(X,VP,R,D,W) 
{
  x = colSums(X[max(VP-R,1):min(VP+R,Xn),max(VP-D,1):min(VP+D,Xn)])
  x = rollsum(x, k=W, align="center", fill=NA)
  return (x)
}

# Create viewpoints list
vp_list = as.numeric(apply(vp_table,1,function(v) { if (v["strand"]=='+') { v["start"] } else { v["end"] } })) %/% U
n_vp = length(vp_list)

# Create target anchors list and initialize output virtual5C file
anchor_list = as.numeric(apply(anchor_table,1,function(v) { if (v["strand"]=='+') { v["start"] } else { v["end"] } })) %/% U
n_anchor = length(anchor_list)
file_v5C = paste0(outdir,'/virtual-5C.csv')
write(file=file_v5C, 
      paste(
        "Source anchor label", 
        "Target anchor label", 
        "Chromosome", 
        "Source anchor position", 
        "Target anchor position", 
        "Anchor distance", 
        "Count", 
        "CPK2B",
        "VP total count",
        "Count sum",
        sep=',' 
      ))

# create distribution table (will be used later for statistical analysis)
distribution=data.frame(distance=seq(min.d.DB,max.d.DB,bin.DB),counts=0,stringsAsFactors = F)

# Process each viewpoint separately
write(paste("Testing",length(vp_list),"viewpoints..."),stderr())
options(scipen=999)                                                       # disable scientific notation
for (k in 1:n_vp) 
{
  write(paste0(chrname," ",round(100*k/n_vp,2),"%"),stderr())
  
  # generate raw virtual 4Cs
  VP = vp_list[k]
  x0 = v4C(X,VP,R,D,W)                                                    # raw counts
  x = round(x0/CPK2B,3)                                                   # normalized counts
  
  # generate virtual 5C data (1)
  vp_offset = VP - max(VP-D,1) + 1                                        # position of viewpoint (VP) in v4C vector
  delta = anchor_list - VP                                                # distances of all anchors from viewpoint
  J = abs(delta)<=D                                                       # indexes of anchors that are within distance D from VP
  anchor_labels = as.character(anchor_table$label[J])                     # labels of these anchors
  anchor_offsets = vp_offset + delta[J]                                   # positions of these anchors in v4C vector
  h = cbind(anchor_offsets - R,anchor_offsets + R)                        # anchors start & end positions in v4C vector (+-R)
  h[h[,1]<0,1]=1
  g = unlist(apply(h,1,function(v) { (v[1]:v[2]) } ))                     # all the positions of the anchor regions (+-R)
  g = g[order(g)]
  
  # get distance
  if(nullRmvAnchors){ a=x[!x %in% g] } else { a=x } 			  # remove target-anchor regions from the null distribution
  d2VP=abs(as.numeric(names(a))-(VP*100))
  d2VP.df=data.frame(d2VP,a)
  names(d2VP.df)=c("distance","counts")
  d2VP.df=d2VP.df[!is.na(d2VP.df$counts),]
  
  df.k=data.frame(distance=seq(min.d.DB,max.d.DB,bin.DB),counts=0,stringsAsFactors = F)
  df.k$start=df.k$distance-r+1
  df.k$end=df.k$distance+r
  
  # get binned count frequency (for statistical analysis)
  i=1
  for (i in 1:nrow(df.k)){ 
    xi=sum(d2VP.df$counts[d2VP.df$distance >= df.k$start[i] & d2VP.df$distance <= df.k$end[i]])
    if (!is.na(xi)){
      df.k$counts[i]=xi
    } else { df.k$counts[i]=0 }
  }
  vp.counts=sum(df.k$counts)                                      # store total viewpoint contacts
  distribution$counts=distribution$counts+df.k$counts             # aggregate distance/counts in master table 
  
  # generate coordinates
  coord_start = as.numeric(names(x))
  coord_end = coord_start + U
  
  if (v4c_bdg){
  	# generate virtual 4C bedgraph files
	if (normalize_bdg) { x_out = cbind(coord_start,coord_end,x) } else { x_out = cbind(coord_start,coord_end,x0) }

  	x_out = x_out[!is.na(x_out[,3]),]
  	x_out = cbind(chrname,x_out)
  	filename = paste(outdir,'/',vp_table$label[k],'-',chrname,'-v4C.bedgraph',sep='') 
  	cat(paste("track type=bedGraph name=",vp_table$label[k],"-",chrname,"\n",sep=""),file=filename)
  	write.table(file=filename,x_out,append=T,quote=F,col.names=F,row.names=F,sep='\t') 
  }

  # generate virtual 5C data (2)
  offsets_sum = as.numeric(apply(h,1,function(v) { sum(x[v[1]:v[2]]) } ))      # sum all the anchor v4c scores in radius  
  K = x[anchor_offsets]>=minvalue                                              # find anchors with enough counts
  K[is.na(K)] = FALSE
  if (sum(K)>0) {
    anchor_data =
      cbind(as.character(vp_table$label[k]),
            anchor_labels,
            chrname,
            coord_start[vp_offset],
            coord_start[anchor_offsets],
            coord_start[anchor_offsets]-coord_start[vp_offset],
            x0[anchor_offsets],
            x[anchor_offsets],
            vp.counts,
            offsets_sum
      )
    anchor_data = anchor_data[K,,drop=FALSE]
    write.table(file=file_v5C, append=TRUE, sep=',', quote=F, row.names=F, col.names=F, anchor_data)
  }
  
}

distribution$chr=chrname
file_distri = paste0(outdir,'/distribution_',chrname,'.tsv')
write.table(x = distribution,file=file_distri, sep='\t', quote=F, row.names=F, col.names=F)

write("Done.",stderr())

quit(save='no')
