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
  make_option(c("--mincount"),default=10, help="minimum count filter for virtual 5C")
)

# process command line arguments
arguments = parse_args(args=commandArgs(trailingOnly=T), OptionParser(usage=usage,option_list=option_list), positional_arguments=c(0,Inf))
opt = arguments$options
inputs = arguments$args
if (length(inputs) != 3) { write("Error: wrong number of inputs! Use --help to see help information", stderr()); quit(save='no') }

# input parameters
outdir = inputs[1]
chrname = inputs[2]
n_reads = as.integer(opt$"nreads")        # number of sequenced read pairs in sample 1
U = as.integer(opt$"unit")                # maximum resolution (bp)
d = as.integer(opt$"maxdist")             # maximum distance from viewpoint
r = as.integer(opt$"radius")              # radius around viewpoint
mincount = as.numeric(opt$"mincount")     # minimum count filter for virtual 5C

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
		 sep=',' 
		))

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
  
  # generate coordinates
  coord_start = as.numeric(names(x))
  coord_end = coord_start + U
  
  # generate virtual 4C bedgraph files
  x_out = cbind(coord_start,coord_end,x)  
  x_out = x_out[!is.na(x_out[,3]),]
  x_out = cbind(chrname,x_out)
  filename = paste(outdir,'/',vp_table$label[k],'-',chrname,'-v4C.bedgraph',sep='') 
  cat(paste("track type=bedGraph name=",vp_table$label[k],"-",chrname,"\n",sep=""),file=filename)
  write.table(file=filename,x_out,append=T,quote=F,col.names=F,row.names=F,sep='\t') 
  
  # generate virtual 5C data
  vp_offset = VP - max(VP-D,1) + 1                                         # position of viewpoint (VP) in v4C vector
  delta = anchor_list - VP                                                 # distances of all anchors from viewpoint
  J = abs(delta)<=D                                                        # indexes of anchors that are within distance D from VP
  anchor_labels = as.character(anchor_table$label[J])                      # labels of these anchors
  anchor_offsets = vp_offset + delta[J]                                    # positions of these anchors in v4C vector
  K = x[anchor_offsets]>=mincount                                          # find anchors with enough counts
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
		    x[anchor_offsets]
		)
    anchor_data = anchor_data[K,,drop=FALSE]
    write.table(file=file_v5C, append=TRUE, sep=',', quote=F, row.names=F, col.names=F, anchor_data)
  }
 
}

write("Done.",stderr())

quit(save='no')

