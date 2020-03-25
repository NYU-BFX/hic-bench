argv = commandArgs(trailingOnly = TRUE)
HKcounts = argv[1L]
TSScounts = argv[2L]
sampleName = argv[3L]

library(ggplot2)

HK.counts=read.csv(HKcounts,sep = "\t",header = F,skip = 1,stringsAsFactors = F)
TSS.counts=read.csv(TSScounts,sep = "\t",header = F,skip = 1,stringsAsFactors = F)
#ActiveMark.counts=read.csv(TSScounts,sep = "\t",header = F,skip = 1, stringsAsFactors = F)

pc1.counts=cbind(HK.counts,TSS.counts[,5])
cols=c("chr","start","end","pc1","HK.counts","TSS.counts")
colnames(pc1.counts)=cols

pc1.counts$group=1
pc1.counts$group[pc1.counts$pc1<0]=-1
pc1.counts$size=pc1.counts$end-pc1.counts$start

summary.df=data.frame(chromosome=rep(unique(pc1.counts$chr),each=3),signGroup=NA,HK.counts=NA,TSS.counts=NA,size.bp=NA,size.per.count.HK=NA,size.per.count.TSS=NA,invertSign=NA,stringsAsFactors = F)
summary.df=summary.df[order(as.numeric(gsub(summary.df$chromosome,pattern="chr",replacement = ""))),]
df.gg = data.frame(ratioMinusPlus=numeric(5*length(unique(summary.df$chromosome))),metric=NA,chr=as.character(rep(unique(pc1.counts$chr),each=5)),stringsAsFactors = F)
suppressWarnings({df.gg=df.gg[order(as.numeric(gsub(df.gg$chr,pattern="chr",replacement = ""))),]})

options(scipen=999)
for (chrom in unique(summary.df$chr)){
  print(chrom)
  pc1.counts.chr=pc1.counts[pc1.counts$chr == chrom,]
  countsBySignHK=aggregate(pc1.counts.chr$HK.counts, list(pc1.counts.chr$group), sum)
  countsBySignTSS=aggregate(pc1.counts.chr$TSS.counts, list(pc1.counts.chr$group), sum)
  sizeBySign=aggregate(pc1.counts.chr$size, list(pc1.counts.chr$group), sum)
  colnames(countsBySignHK)=c("signGroup","HK.counts")
  colnames(countsBySignTSS)=c("signGroup","TSS.counts")
  colnames(sizeBySign)=c("signGroup","size.bp")
  info.df=merge(countsBySignHK,countsBySignTSS)
  info.df=merge(info.df,sizeBySign)
  info.df$size.per.count.HK=info.df$size.bp/info.df$HK.counts 
  info.df$size.per.count.TSS=info.df$size.bp/info.df$TSS.counts 
  ratioMinusPlus = log(info.df[2,]/info.df[1,]) #ratio of 1 / -1 
  info.df=rbind(info.df,ratioMinusPlus)
  #info.df=rbind(info.df,info.df[3,1:ncol(info.df)]-1 ) 
  info.df[3,1]=c("log(ratio[1:-1])")
  df.gg$ratioMinusPlus[df.gg$chr == chrom] = unlist(ratioMinusPlus[2:ncol(ratioMinusPlus)])
  df.gg$metric[df.gg$chr == chrom] = colnames(ratioMinusPlus[2:ncol(ratioMinusPlus)])
  
  if(info.df$HK.counts[1] & info.df$TSS.counts[1] > 0){
    
    if (info.df$size.per.count.HK[3] > 1 & info.df$size.per.count.TSS[3] > 1){ isA = -1 } 
    if (info.df$size.per.count.HK[3] < 1 & info.df$size.per.count.TSS[3] < 1){ isA = 1} 
    if (info.df$size.per.count.HK[3] < 1 & info.df$size.per.count.TSS[3] > 1 | info.df$size.per.count.HK[3] > 1 & info.df$size.per.count.TSS[3] < 1){ 
      x=as.data.frame(abs(info.df[3,5:6])==max(abs(info.df[3,5:6]),na.rm = T))
      isMaxRatio=colnames(x)[grep("TRUE",x,value = F)]
      if(info.df[3,isMaxRatio] < 0){isA = 1} else { isA = -1}
    }
  } 
  
  summary.df[summary.df$chr == chrom,2:(ncol(summary.df)-1)]=info.df
  
  ## Fix pc1 sign if required ##
  print(paste0("Compartment A is: ",isA))
  if (isA < 0){ 
    pc1.counts$pc1[pc1.counts$chr == chrom]=pc1.counts$pc1[pc1.counts$chr == chrom]*-1 
    print("pc1 sign inverted")
    summary.df$invertSign[summary.df$chromosome==chrom]=T
  } else {
    print("Keep pc1 sign")
    summary.df$invertSign[summary.df$chromosome==chrom]=F
  }
  
}

summary.df[,3:(ncol(summary.df)-1)]=apply(summary.df[,3:(ncol(summary.df)-1)],2,function(x){round(x,digits=3)})
summary.df$sampleName=sampleName
#write.table(pc1.counts[,1:4],"pca_HKgenesFix.PC1.bedGraph",quote = F,col.names = F,row.names = F,sep="\t")
write.table(summary.df,"pc1_metrics_summary.txt",quote = F,col.names = T,row.names = F,sep="\t")

##### Downstream Analysis ####
df.gg$chr=as.factor(df.gg$chr)
df.gg=df.gg[df.gg$chr != "chrY",]
df.gg$chr=factor(df.gg$chr,levels=unique(df.gg$chr))

##1 (raw)
#pdf('metrics_AB.ratio_all.pdf',width = 12,useDingbats=FALSE)
#ggplot(df.gg,aes(chr,ratioMinusPlus,color=metric))+
#  geom_point(size=3)+
#  ylab("log ratio [A:B]")+
#  xlab("")+
#  geom_hline(yintercept = 0,color="darkred", linetype="dashed")
#+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
#dev.off()

##2 (turn ratio to a score by flipping the sign of negatively correlated metrics)
df.gg2=df.gg
df.gg2$ratioMinusPlus[df.gg2$metric %in% c("size.per.count.HK","size.per.count.TSS", "size.bp")]=df.gg2$ratioMinusPlus[df.gg2$metric %in% c("size.per.count.HK","size.per.count.TSS", "size.bp")]*-1

#pdf('metrics_AB.ratio_sameDirection.pdf',width = 12,useDingbats=FALSE)
#ggplot(df.gg2,aes(chr,ratioMinusPlus,color=metric))+
#  geom_point(size=3)+
#  ylab("log ratio [A:B]")+
#  geom_hline(yintercept = 0,color="darkred", linetype="dashed")
#xlab("")
#dev.off()

##3
df.gg3=df.gg2
df.gg3=df.gg3[df.gg3$metric %in% c("size.per.count.HK","size.per.count.TSS","size.bp"),]

pdf('metrics_AB.logratio.pdf',width = 12,useDingbats=FALSE)
ggplot(df.gg3,aes(chr,ratioMinusPlus,color=metric))+
  geom_point(size=3)+
  ylab("log ratio [A:B]")+
  xlab("")+
  geom_hline(yintercept = 0,color="darkred", linetype="dashed")
dev.off()
