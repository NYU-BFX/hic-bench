
v4C = function(X,r,R,W) 
{
  x = colSums(X[(r-R):(r+R),])
  x = rollsum(x, k=W, align="center")
  return (x)
}

X = readMM('results/matrix-sparse.by_sample.maxd_5Mb/filter.by_sample.mapq_20/align.by_sample.bowtie2/BALL-9_R-Arima-rep1/matrix.chr15.mtx')

U = 100
R = 5000 %/% U
W = 20000 %/% U
D = 250000 %/% U

vp = 63434260 %/% U
I = (vp-D):(vp+D)
X1 = X[I,I]

n = (nrow(X1)-2*R)
Y1 = matrix(0,n,n)
for (k in 1:n) 
{
  v = v4C(X1,R+k,R,W)
  j = (n-length(v)) %/% 2
  Y1[k,j:(length(v)+j-1)] = v
}

Y2 = Y1
z = abs(row(Y2)-col(Y2))<=W
Y2[z] = 0

png('a.png'); image(Y2); dev.off()



library(raster)

## Convert matrix to a raster object
r <- raster(Y2)
extent(r) <- extent(c(0, nrow(Y2), 0, nrow(Y2)) + 0.5)

## Find the maximum value within the 9-cell neighborhood of each cell
f <- function(X) max(X, na.rm=TRUE)
ww <- matrix(1, nrow=10, ncol=10) ## Weight matrix for cells in moving window
localmax <- focal(r, fun=f, w=ww, pad=TRUE, padValue=NA)

## Does each cell have the maximum value in its neighborhood?
r2 <- r==localmax

png('a.png'); image(Y2); dev.off()

png('b.png'); image(r2); dev.off()

