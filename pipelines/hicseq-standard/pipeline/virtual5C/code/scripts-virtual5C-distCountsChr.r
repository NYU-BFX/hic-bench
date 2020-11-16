## generates a distribution table of counts per distance ##
## per chromosome on a random set of genomic regions ##
## example: Rscript distCounts.r ./distribution/ chr16 100000000 100 2000000 2500 1000 hg19_genome.bed matrix.chr16.mtx 

# functions
getRawCounts = function(X,VP,R,D){        # get raw counts from matrix sparse X in +-D window from VP site at 2*R resolution
  x = colSums(X[max(VP-R,1):min(VP+R,Xn),max(VP-D,1):min(VP+D,Xn)])
  return(x)
}

# load libraries
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(zoo))

# set parameters
argv = commandArgs(trailingOnly = T)
outdir = argv[1L]
chrname = argv[2L]
n_reads = as.numeric(argv[3L])                        # number of intrachromosomal valid reads
U = as.numeric(argv[4L])                              # maximum matrix sparse resolution (bp)
d = as.numeric(argv[5L])                              # maximum distance from viewpoint (vp = random site)
r = as.numeric(argv[6L])                              # radius (resolution * 0.5)
tsize = as.numeric(argv[7L])                          # number of random sites per chromosome used for the distribution calculation
chrdata = argv[8L]                         # chromosome length data file
mat1 = argv[9L]                            # hic chr-matrix 

# adjust by unit 
R = r %/% U           # number of bins to extend around each side of peak center (100bp) to anchor the vp  (i.e=25)
D = d %/% U           # total number of bins to extend each side of peak center (100bp) to test the vp (i.e=10000)
W = 2*R + 1           # total bins covering each vp and each anchor (resolution: i.e: 50)

# adjust maximum distance (for consistency)
D = D + R*2

# load chromosome length data
chr_data = read.table(chrdata,stringsAsFactors = F)
names(chr_data) = c("chr","start","end","chr.id")

# get random sites
chrl = chr_data$end[chr_data$chr==chrname]
chr_range = (d*2):(chrl-(d*2))
vp_list  =as.numeric(sample(chr_range,tsize)) %/% U
n_vp = length(vp_list)
  
# load matrix
X = readMM(mat1)                    # raw counts per 100bp bins around the peak center (+-1Mb)
X = X+t(X); diag(X) = diag(X)/2     # make symmetric
colnames(X) = seq(0,ncol(X)-1)*U    # set column labels by coordinate start
Xn = nrow(X)
  
# compute and store contact vs distance distribution for each random site
l = (D*2)+1
distCounts = matrix(nrow = l*tsize, ncol=2)
rs = 1    # row start index
re = l    # row end index
  
for (k in 1:n_vp){
  write(paste0(chrname," ",round(100*k/n_vp,2),"%"),stderr())
  VP = vp_list[k]
  raw_counts = getRawCounts(X,VP,R,D)                   # i.e. raw counts per 100bp bins around the peak center (+-1Mb)
  d2VP=abs(as.numeric(names(raw_counts))-(VP*100))
  distCounts[rs:re,]=cbind(d2VP,raw_counts)
  rs=rs+l
  re=re+l
}

distCounts = as.data.frame(distCounts)
names(distCounts) = c("distance","counts")
distCounts$chr = chrname

outname = paste0(outdir,'/distCounts_',chrname,'.tsv')
write.table(x = distCounts,file=outname, sep='\t', quote=F, row.names=F, col.names=F)

write("Done.",stderr())
