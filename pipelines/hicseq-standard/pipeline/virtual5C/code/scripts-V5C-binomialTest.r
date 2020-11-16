### FUNCTIONS ###
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
# set parameters
argv = commandArgs(trailingOnly = T)
spline = argv[1L]
v5c_file = argv[2L] 
outdir = argv[3L]
sp_res = as.numeric(argv[4L])
r = as.numeric(argv[5L])                              # radius around viewpoint
min.d.DB = as.numeric(argv[6L])			                  # minimum distance computed in the null distribution
max.d.DB = as.numeric(argv[7L])			                  # maximum distance computed in the null distribution
bin.DB = as.numeric(argv[8L])				                  # bin size used to compute the null distribution

# spline = "/Users/javrodher/Desktop/spline_corrected.csv"
# r=2500
# min.d.DB=20000
# max.d.DB=1000000
# bin.DB=5000
# sp_res = 100
# outdir="/Users/javrodher/Desktop/"
# v5c_file = "/Users/javrodher/Work/RStudio-PRJs/HICBENCH-DEV/V5C_backgroundMod/virtual-5C.csv"				            # V5C data

# load corrected spline
corrected_spline = read.csv(spline,stringsAsFactors = F)
N=corrected_spline$N[1]

# perform binomial test
v5c_data=read.csv(v5c_file,stringsAsFactors = F)
v5c_data=v5c_data[v5c_data$Anchor.distance > min.d.DB,]
v5c_data=v5c_data[!is.na(v5c_data$Count),]

v5c_data = binomTest(v5c_data,N,corrected_spline) # binomial test computation

ttl = nrow(v5c_data)
tsl = sum(v5c_data$fdr< 0.01)
tslp = round(tsl*100/ttl)

write(paste0("Total tested loops = ",ttl),stderr())
write(paste0("Total significant loops = ",tsl," (",tslp,"%)"),stderr())
options(scipen = 900)

# save
if (nrow(v5c_data) > 0){
  v5c_data[,9:12]=apply(v5c_data[,9:12],2,function(x) formatC(x,format = "e", digits = 2))
  write.table(v5c_data,paste0(outdir,"/virtual-5C_btest.csv"), sep=',', quote=F, row.names=F, col.names=T) 
}
