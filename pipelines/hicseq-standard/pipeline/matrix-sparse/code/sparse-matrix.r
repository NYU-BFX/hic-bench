#!/usr/bin/Rscript
#$ -S /usr/bin/Rscript

args = commandArgs(trailingOnly=TRUE)

chr_matrix = args[1]         # e.g. TALL/matrix.chr8.mtx

library(Matrix)
library(zoo)

vp = 128747600     # viewpoint coordinate (bp)
r1 = 0             # radius around viewpoint (bp)
r2 = 10000         # radius around viewpoint (bp)
d = 2500000        # maximum distance from viewpoint (bp)
w = 100000         # smoothing window (bp)

# adjust by unit 
U = 100            # bp
VP = vp %/% U
R1 = r1 %/% U
R2 = r2 %/% U
D = d %/% U
W = w %/% U

# load matrix
X <- readMM(chr_matrix)
X = X+t(X)
diag(X) = diag(X)/2
print(object.size(X), units = "auto", standard = "SI")

# virtual 4C
s = rollsum(colSums(X[(VP-R1):(VP+R2),(VP-D):(VP+D)]),k=W)

# write output
write.table(file='x.csv',s,col.names=F,row.names=F,quote=F)

quit(save='no')

