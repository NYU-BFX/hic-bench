### FUNCTIONS ###
getDistribution = function(x,M,N){
  # average counts by distance bin + spline + antitonic correction
  X=x[x$distance>=min.d.DB-r & x$distance<=max.d.DB+r,]
  dbins=seq(min.d.DB,max.d.DB,bin.DB)
  dtab.random=data.frame(distance=dbins,counts.sum=0,counts.avg=0,davg=NA,stringsAsFactors = F)
  dtab.random$start=dtab.random$distance-r+1
  dtab.random$end=dtab.random$distance+r
  
  for (i in 1:nrow(dtab.random)){
    print(i)
    xi=X[X$distance >= dtab.random$start[i] & X$distance <= dtab.random$end[i],]
    xi=xi[!is.na(xi$counts),]
    davg=mean(xi$distance[xi$distance >= dtab.random$start[i] & xi$distance <= dtab.random$end[i]])
    csum=sum(xi$counts[xi$distance >= dtab.random$start[i] & xi$distance <= dtab.random$end[i]])
    cavg=mean(xi$counts[xi$distance >= dtab.random$start[i] & xi$distance <= dtab.random$end[i]])
    dtab.random$counts.sum[i]=csum
    dtab.random$counts.avg[i]=cavg
    dtab.random$davg[i]=csum=davg
  }
  
  dtab.random$hi=dtab.random$counts.sum/M         # avg contacts per unique locus-pair (expected contacts)
  dtab.random$pe=dtab.random$hi/N                 # hi / N  # expected probability
  
  smsp=spline(x=dtab.random$davg,                 # fit spline
              y=dtab.random$pe,
              n=nrow(dtab.random)*(bin.DB/sp_res))
  
  stats.spl=data.frame(distance=smsp$x,           # make spline data frame
                       pe=smsp$y,
                       stringsAsFactors = F)
  
  mr = monoreg(stats.spl$distance,stats.spl$pe, type="a")     # correct spline so it is monotonic descendant
  stats.spl_antitn=data.frame(distance=mr$x,
                              pe=mr$yf,
                              stringsAsFactors = F)
  
  return(stats.spl_antitn)
}

# compute binomial test
binomTest = function(x,N,nullData){
  options(scipen = 1)
  trials=N
  x$pvalue=NA
  x$p.expected=NA
  x$p.observed=NA
  
  for (i in 1:nrow(x)){
    write(paste0(round(100*i/nrow(x),2),"%"),stderr())
    d=abs(x$Anchor.distance[i])
    successes=round(x$Count[i])
    temp=nullData
    temp$d=abs(nullData$distance-d)
    p.null=mean(temp$pe[temp$d==min(temp$d)[1]])
    
    if(trials>0) { 
      test=binom.test(successes, trials, p.null,alternative="greater")
      x$pvalue[i]=test$p.value
    } else { x$pvalue[i]=1 }
    x$p.expected[i]=p.null
    x$p.observed[i]=successes/trials
  }
  x$fdr=p.adjust(x$pvalue,n=nrow(x),method="fdr")
  return(x)
}

### START ###
# load libraries
suppressPackageStartupMessages(library(fdrtool))
suppressPackageStartupMessages(library(ggplot2))

# set parameters
argv = commandArgs(trailingOnly = T)
distCounts = argv[1L]
outdir = argv[2L]
sp_res = as.numeric(argv[3L])
r = as.numeric(argv[4L])                              # radius around viewpoint
min.d.DB = as.numeric(argv[5L])			                  # minimum distance computed in the null distribution
max.d.DB = as.numeric(argv[6L])			                  # maximum distance computed in the null distribution
bin.DB = as.numeric(argv[7L])				                  # bin size used to compute the null distribution
tsize = as.numeric(argv[8L])	
nchroms = as.numeric(argv[9L])	

# print(distCounts)
# print(outdir)
# print(sp_res)
# print(r)
# print(min.d.DB)
# print(max.d.DB)
# print(bin.DB) 
print(tsize)
print(nchroms)

# distCounts = "/Users/javrodher/Work/RStudio-PRJs/HICBENCH-DEV/data/distribution.tsv"
# r=2500
# min.d.DB=20000
# max.d.DB=1000000
# bin.DB=5000
# sp_res = 100
# outdir="/Users/javrodher/Desktop/"
# tsize=50
# nchroms=20

# load distribution data
data1 = read.delim(distCounts,header = F)
names(data1)=c("distance","counts","chr")
#data1=data1[sample(1:nrow(data1),1000000),]  # for quick testing

# compute raw null distribution
M = tsize*2*nchroms         # total number of unique locus-pairs tested
N = sum(data1$counts)       # total number of counts in locus-pairs tested

spline_raw = getDistribution(data1,M,N)       # average counts by distance bin + spline + antitonic correction
spline_raw$N = N                    # save N parameters here to be used in the final binomial test

spline_raw_sparse = spline_raw
spline_raw_sparse$pe = spline_raw_sparse$pe/(bin.DB/sp_res)

# compute corrected null distribution
outlThr = 1/M           # outlier threshold (NOT WORKING)
data1$id = 1:nrow(data1)

data2 = data1[data1$counts>0,]
names(data2)[1:2] = c("Anchor.distance","Count")

# make sure there arent any probability value outside of 0-1 range (spline doesn't care about this)
spline_raw_sparse$pe[spline_raw_sparse$pe < 0 ] = 0
spline_raw_sparse$pe[spline_raw_sparse$pe > 1 ] = 1

options(scipen = 1)
data2 = binomTest(data2,N,spline_raw_sparse)  # binomial test computation

#hist(data2$Anchor.distance[data2$pvalue <= outlThr],breaks=1000)
#hist(data2$Anchor.distance[data2$fdr <= outlThr],breaks=1000)
#plot(data2$Anchor.distance,data2$pvalue)

tol = sum(data2$pvalue <= 0.0001)
write(paste0("Total outliers = ",tol),stderr())

data2 = data2[data2$pvalue > 0.0001,]            
data_clean = data1[!data1$id %in% data2$id,]
Ncorrected = sum(data_clean$counts)        
spline_corrected = getDistribution(data_clean,M,Ncorrected)          
spline_corrected$N = Ncorrected                    # save N parameters here to be used in the final binomial test

# save data
write.csv(spline_corrected,paste0(outdir,"/spline_corrected.csv"),row.names = F)

# plot raw and corrected splines
pdf(paste0(outdir,"/splines.pdf"), width=6 ,height=6) 
ggplot()+
  geom_line(data=spline_corrected,aes(x=distance,y=pe),size=1,color="blue")+
  geom_line(data=spline_raw,aes(x=distance,y=pe),size=1,color="red")
dev.off()

