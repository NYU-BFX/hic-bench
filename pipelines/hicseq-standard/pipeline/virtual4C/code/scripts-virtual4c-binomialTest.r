argv = commandArgs(trailingOnly = T)
distribution_file = argv[1L]
v5c_file = argv[2L]
sparse_res = argv[3L]
fdr.cut = argv[4L]
p.cut = argv[5L]
outdir = argv[6L]

library(ggplot2)
# parameters
exclude_chr = c("chrX","chrY")
options(scipen = 900)
cols = c("chr","start","end")

# load data
null_data=read.table(distribution_file,stringsAsFactors = F)
v5c_data=read.csv(v5c_file,stringsAsFactors = F)
head(v5c_data)
names(null_data)=c("distance","counts.sum","chr")

# frequency -> probability
null_data=null_data[!null_data$chr %in% exclude_chr,]
null_data$p=null_data$counts.sum/sum(null_data$counts.sum)

# unify chromosome data by distance
stats=data.frame(distance=unique(null_data$distance),p=NA,stringsAsFactors = F)
for (d in unique(stats$distance)){ stats$p[stats$distance==d]=sum(null_data$p[null_data$distance==d]) }
stats=stats[-nrow(stats),]

# spline
stats.spl=data.frame(spline(stats$distance,stats$p,n = nrow(stats)))
names(stats.spl)=c("distance","p")
stats.spl$p=stats.spl$p/sum(stats.spl$p)
stats.spl$distance=round(stats.spl$d)

# binomial test
v5c_data=v5c_data[abs(v5c_data$Anchor.distance) > stats.spl$distance[1] & abs(v5c_data$Anchor.distance) < stats.spl$distance[nrow(stats.spl)],]
v5c_data=v5c_data[!is.na(v5c_data$Count.sum),]
v5c_data$pvalue=NA
v5c_data$p.expected=NA

for (i in 1:nrow(v5c_data)){
  trials=round(v5c_data$VP.total.count[i])
  d=abs(v5c_data$Anchor.distance[i])
  successes=round(v5c_data$Count.sum[i])
  temp=stats.spl
  temp$d=abs(stats.spl$distance-d)
  p.null=mean(temp$p[temp$d==min(temp$d)[1]])
  if(trials>0) { 
    test=binom.test(successes, trials, p.null,alternative = "greater")
    v5c_data$pvalue[i]=test$p.value
  } else { v5c_data$pvalue[i]=1 }
  v5c_data$p.expected[i]=p.null
}
v5c_data$fdr=p.adjust(v5c_data$pvalue,n = nrow(v5c_data),method = "fdr")
v5c_data$p=v5c_data$Count.sum/v5c_data$VP.total.count

# plots
# distribution spline
pdf(paste0(outdir,"/spline_null_distribution.pdf"), width=6 ,height=6) 
ggplot(data=stats.spl,aes(x=distance,y=p))+
  geom_point(color="blue")+
  ylab("expected contact probability")+
  xlab("genomic distance (bp)")
dev.off()

# reorganize v5c table
names(v5c_data)[14]="p.observed"
v5c_data=v5c_data[,c(1:11,13,12,14)]

# save
write.table(v5c_data,paste0(outdir,"/virtual-5C.csv"), sep=',', quote=F, row.names=F, col.names=T)
