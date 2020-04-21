argv = commandArgs(trailingOnly = TRUE)
inpdir = argv[1L]
n_chr = argv[2L]
outdir = "./QC_plots/"

library(ggplot2)

setwd(inpdir)
chroms=paste0("chr",1:as.numeric(n_chr))
chroms=c(chroms,"chrX")
infiles=list.files(pattern="filtered")
infiles=infiles[grep(".tsv",infiles)]
infiles=infiles[grep("raw",infiles)]
infiles=infiles[grep(".gz",infiles,invert = T)]
labels=gsub(infiles,pattern = ".tsv",replacement = "")
options(scipen=10000)

for(file.num in 1:length(infiles)){
  loops.df=read.table(infiles[file.num],header = T,stringsAsFactors = F)
  loops.df=loops.df[loops.df$chr1 != "chrY",]
  loops.df$distance=loops.df$fragmentMid2-loops.df$fragmentMid1
  df=loops.df
  df$q.value=-log2(df$q.value)
  df$distance=abs(df$distance/1000000)
  df$chr=as.numeric(gsub(df$chr1,pattern = "chr",replacement = ""))
  df=df[order(df$chr),]
  df$chr[is.na(df$chr)]="X"
  df$chr1=as.factor(df$chr1)
  df$chr1=factor(df$chr1,levels=chroms)
  missing.chr=chroms[!chroms %in% df$chr1]
  for (i in missing.chr){ df=rbind(df,c(i,0,0,0,0,0,0,0,0,0,0))}
  df$q.value=as.numeric(df$q.value)
  df$distance=as.numeric(df$distance)
  
  pos=data.frame(x=0,y=c(0.001,0.01,0.05))
  pos$label=pos$y
  pos$y=-log2(pos$y)
  
  pos.chr=pos
  pos.chr$x=15
  pos.chr$label=c("***","**","*")
  
  pos.histo=pos.chr
  pos.histo$x=0
  bias.df=data.frame(bias=as.numeric(c(df$bias1,df$bias2)),chr1=rep(df$chr1,2),stringsAsFactors = F)
  
  ## Generate QC report ##
  pdf(paste0(outdir,labels[file.num],"_QC.pdf"),width=12, height=8)
  # bias histogram
  print(ggplot(data=bias.df,aes(x= bias))+
          geom_histogram()+
          xlim(-2,2)+
          geom_vline(xintercept = 1,linetype="dotted",color="red")+
          geom_vline(xintercept = -1,linetype="dotted",color="red")+
          ggtitle("loops bias-score distribution")+
          theme(plot.title = element_text(color = "black",hjust = 0.5))+
          ylab("loop-anchor frequency")+
          xlab("bias score"))
  
  # bias histogram by chromosome
  print(ggplot(data=bias.df,aes(x= bias))+
          geom_histogram()+
          xlim(-2,2)+
          geom_vline(xintercept = 1,linetype="dotted",color="red")+
          geom_vline(xintercept = -1,linetype="dotted",color="red")+
          ggtitle("loops bias-score distribution by chromosome")+
          theme(plot.title = element_text(color = "black",hjust = 0.5))+
          facet_wrap(chr1~.,scales = "fixed",nrow=3)+
          ylab("loop-anchor frequency")+
          xlab("bias score"))
  
  #ALL
  print(ggplot(data=df,aes(x = distance,y = q.value))+
          geom_smooth()+
          xlim(c(0,20))+
          geom_hline(yintercept = -log2(0.05),color="red",linetype = "dotted")+
          geom_hline(yintercept = -log2(0.01),color="red",linetype = "dotted")+
          geom_hline(yintercept = -log2(0.001),color="red",linetype = "dotted")+
          ggtitle("loops qvalue-distance smoothed distribution")+
          theme(plot.title = element_text(color = "black",hjust = 0.5))+
          ylab("-log2(qvalue)")+
          xlab("distance (Mb)")+
          geom_text(data=pos,aes(x=x,y=y,label=label),size=3,vjust=-1))
  
  #CHROMS
  print(ggplot(data=df,aes(x = distance,y = q.value))+
          geom_smooth()+
          xlim(c(0,20))+
          facet_wrap(chr1~., scales= "fixed",nrow=3)+
          geom_hline(yintercept = -log2(0.05),color="red",linetype = "dotted")+
          geom_hline(yintercept = -log2(0.01),color="red",linetype = "dotted")+
          geom_hline(yintercept = -log2(0.001),color="red",linetype = "dotted")+
          ggtitle("loops qvalue-distance smoothed distribution by chromosome")+
          theme(plot.title = element_text(color = "black",hjust = 0.5))+
          ylab("-log2(qvalue)")+
          xlab("distance (Mb)")+
          geom_text(data=pos.chr,aes(x=x,y=y,label=label),size=2,vjust=0,hjust=0.5))
  
  # qvalue histogram
  print(ggplot(data=df,aes(x= q.value))+
          geom_histogram()+
          xlim(c(-1,20))+
          geom_vline(xintercept = -log2(0.05),color="red",linetype = "dotted")+
          geom_vline(xintercept = -log2(0.01),color="red",linetype = "dotted")+
          geom_vline(xintercept = -log2(0.001),color="red",linetype = "dotted")+
          ggtitle("loops qvalue-distance distribution")+
          theme(plot.title = element_text(color = "black",hjust = 0.5))+
          ylab("loop frequency")+
          xlab("-log2(qvalue)")+
          geom_text(data=pos.chr,aes(x=y,y=x,label=label),size=3,vjust=2))
  
  # qvalue histogram by chromosome
  print(ggplot(data=df,aes(x= q.value))+
          geom_histogram()+
          xlim(c(-1,20))+
          facet_wrap(chr1~., scales= "fixed",nrow=3)+
          geom_vline(xintercept = -log2(0.05),color="red",linetype = "dotted")+
          geom_vline(xintercept = -log2(0.01),color="red",linetype = "dotted")+
          geom_vline(xintercept = -log2(0.001),color="red",linetype = "dotted")+
          ggtitle("loops qvalue-distance distribution by chromosome")+
          theme(plot.title = element_text(color = "black",hjust = 0.5))+
          ylab("loop frequency")+
          xlab("-log2(qvalue)")+
          geom_text(data=pos.chr,aes(x=y,y=x,label=label),size=2,vjust=1.2))
  
  # distance histogram
  print(ggplot(data=df,aes(x= distance))+
          geom_histogram()+
          xlim(c(0,10))+
          ggtitle("loops distance distribution")+
          theme(plot.title = element_text(color = "black",hjust = 0.5))+
          ylab("loop frequency")+
          xlab("distance (Mb)"))  
  
  # distance histogram by chromosome
  print(ggplot(data=df,aes(x= distance))+
          geom_histogram()+
          facet_wrap(chr1~., scales= "fixed",nrow=3)+
          xlim(c(0,10))+
          ggtitle("loops distance distribution by chromosome")+
          theme(plot.title = element_text(color = "black",hjust = 0.5))+
          ylab("loop frequency")+
          xlab("distance (Mb)")+
          theme(axis.text.x = element_text(angle = 90)))
  dev.off()
}
