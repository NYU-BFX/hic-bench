argv = commandArgs(trailingOnly = TRUE)
inpdir = argv[1L]
APA_res = as.numeric(argv[2L])
C1 = argv[3L]
C2 = argv[4L]
URm = as.numeric(argv[5L])
outdir = "in-house-analysis"

#inpdir="/Users/javrodher/Work/RStudio-PRJs/leukemia-cell-line-DR/data/loops-diff-may21st/APA/quantiles/"
#APA_res=10000
#C1="CUTLL1_DMSO_A"
#C2="CUTLL1_THZ1"

library(ggplot2)
library(gridExtra)
library(pheatmap)
library(RColorBrewer)
library(grid)
library(ggplotify)

setwd(inpdir)
dir.create(outdir)
options(scipen=10000)
pdf(NULL)

# Functions #
rotate=function(x) {t(apply(x, 2, rev))}

ggHmapAPA=function(mat.path,title,URm){
  if(length(mat.path)==1){
    x=read.csv(mat.path,header = F)
    x$V1=gsub(pattern = "[",replacement = "",x = x$V1,fixed = T)
    x[,ncol(x)]=gsub(pattern = "]",replacement = "",x = x[,ncol(x)],fixed = T)
    x=rotate(x)
    x=apply(x,2,as.numeric)
    UR=mean(x[1:6,(ncol(x)-6):ncol(x)])*URm
    x[x > UR]=UR
    min=min(x)
    n=nrow(x)
    df.hmap=data.frame(row=rep(1:n,each=n),col=rep(1:n,times=n),values=NA,stringsAsFactors = F)
    for(i in 1:nrow(df.hmap)){df.hmap$values[i]=x[df.hmap$row[i],df.hmap$col[i]]}
  } else{ 
    df.hmap=data.frame(row=1,col=1,values=0,stringsAsFactors = F)
    UR=0
  }
  
  print(ggplot(df.hmap, aes(row, col, fill= values)) + 
          geom_tile(show.legend = F,color=NA,size=0)+
          xlab(paste0("bins (",APA_res/1000," kb)"))+
          ylab(paste0("bins (",APA_res/1000," kb)"))+
          ggtitle(title)+
          scale_fill_gradient(low="white",high="red",limits=c(0,UR))+
          theme(plot.title = element_text(hjust = 0.5,size=8),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                panel.background = element_blank(), axis.line = element_blank()))
}  

prepMetrics=function(measures.file,q,s){
  if(length(measures.file)==1){
    x=read.table(measures.file,stringsAsFactors = F)
    m=as.data.frame(t(x[2]))
    names(m)=x$V1
  } else { m=data.frame(P2M=NA,P2UL=NA,P2UR=NA,P2LL=NA,P2LR=NA,ZscoreLL=NA,stringsAsFactors = F) }
  m$q=q
  m$s=s
  return(m)
}

#m=md.all.merge
#metric="P2LL"
#qtab=qcpmtab
#xlab=""
plotMetrics=function(m,metric,xlab,qtab){
  t=read.csv(qtab,stringsAsFactors = F)
  t[,1]=round(t[,1],digits = 3)
  print(ggplot(m,aes(x=q,y=!!ensym(metric),group=group,color=group))+
          geom_point()+
          geom_line()+
          xlab(xlab)+
          geom_hline(yintercept = log2(1),linetype=3,color="red")+
          scale_x_continuous(breaks=c(1:10),labels=paste0("q",1:10),sec.axis = dup_axis(labels=c(t[,1],paste0(">",t[9,1])))))
}

##### RUN #####
all.files=list.files(pattern=".txt",recursive = T,full.names = F,include.dirs = F)
tables=list.files(pattern=".csv",recursive = T,full.names = F,include.dirs = F)

methods=c("APA","rankAPA","centerNormedAPA","normedAPA")
#methods=c("APA")
scores=c("P2M","P2UL","P2UR","P2LL","P2LR","ZscoreLL")
#scores=c("P2LL")

for (method in methods){
  print(method)
  mth=paste0("^",method,".txt")
  
  for (score in scores){
    print(score)
    mats=all.files
    
    qcpmtab=tables[grep("qcpm",tables)]
    qdisttab=tables[grep("qdist",tables)]
    
    #APA plots
    mats=list.files(pattern=mth,recursive = T)
    mats.c1=mats[grep(C1,mats)]
    mats.c2=mats[grep(C2,mats)]
    
    mats.dist.c1=mats.c1[grep("qdist",mats.c1)][c(1,3:10,2)]
    mats.dist.c2=mats.c2[grep("qdist",mats.c2)][c(1,3:10,2)]
    
    mats.q=mats[!mats %in% c(mats.dist.c1,mats.dist.c2)]
    mats.q.c1=mats.q[grep(C1,mats.q)]
    mats.q.c2=mats.q[grep(C2,mats.q)]
    
    mats.q.c1=mats.q.c1[grep("cpm",mats.q.c1)][c(1,3:10,2)]
    mats.q.c2=mats.q.c2[grep("cpm",mats.q.c2)][c(1,3:10,2)]
    
    qdist1.c1=ggHmapAPA(mats.dist.c1[1],"1st quantile",URm)
    qdist2.c1=ggHmapAPA(mats.dist.c1[2],"2nd quantile",URm)
    qdist3.c1=ggHmapAPA(mats.dist.c1[3],"3rd quantile",URm)
    qdist4.c1=ggHmapAPA(mats.dist.c1[4],"4rd quantile",URm)
    qdist5.c1=ggHmapAPA(mats.dist.c1[5],"5th quantile",URm)
    qdist6.c1=ggHmapAPA(mats.dist.c1[6],"6th quantile",URm)
    qdist7.c1=ggHmapAPA(mats.dist.c1[7],"7th quantile",URm)
    qdist8.c1=ggHmapAPA(mats.dist.c1[8],"8th quantile",URm)
    qdist9.c1=ggHmapAPA(mats.dist.c1[9],"9th quantile",URm)
    qdist10.c1=ggHmapAPA(mats.dist.c1[10],"10th quantile",URm)
    
    qdist1.c2=ggHmapAPA(mats.dist.c2[1],"1st quantile",URm)
    qdist2.c2=ggHmapAPA(mats.dist.c2[2],"2nd quantile",URm)
    qdist3.c2=ggHmapAPA(mats.dist.c2[3],"3rd quantile",URm)
    qdist4.c2=ggHmapAPA(mats.dist.c2[4],"4rd quantile",URm)
    qdist5.c2=ggHmapAPA(mats.dist.c2[5],"5th quantile",URm)
    qdist6.c2=ggHmapAPA(mats.dist.c2[6],"6th quantile",URm)
    qdist7.c2=ggHmapAPA(mats.dist.c2[7],"7th quantile",URm)
    qdist8.c2=ggHmapAPA(mats.dist.c2[8],"8th quantile",URm)
    qdist9.c2=ggHmapAPA(mats.dist.c2[9],"9th quantile",URm)
    qdist10.c2=ggHmapAPA(mats.dist.c2[10],"10th quantile",URm)
    
    cpm1.c1=ggHmapAPA(mats.q.c1[1],"1st quantile",URm)
    cpm2.c1=ggHmapAPA(mats.q.c1[2],"2nd quantile",URm)
    cpm3.c1=ggHmapAPA(mats.q.c1[3],"3rd quantile",URm)
    cpm4.c1=ggHmapAPA(mats.q.c1[4],"4th quantile",URm)
    cpm5.c1=ggHmapAPA(mats.q.c1[5],"5th quantile",URm)
    cpm6.c1=ggHmapAPA(mats.q.c1[6],"6th quantile",URm)
    cpm7.c1=ggHmapAPA(mats.q.c1[7],"7th quantile",URm)
    cpm8.c1=ggHmapAPA(mats.q.c1[8],"8th quantile",URm)
    cpm9.c1=ggHmapAPA(mats.q.c1[9],"9th quantile",URm)
    cpm10.c1=ggHmapAPA(mats.q.c1[10],"10th quantile",URm)
    
    cpm1.c2=ggHmapAPA(mats.q.c2[1],"1st quantile",URm)
    cpm2.c2=ggHmapAPA(mats.q.c2[2],"2nd quantile",URm)
    cpm3.c2=ggHmapAPA(mats.q.c2[3],"3rd quantile",URm)
    cpm4.c2=ggHmapAPA(mats.q.c2[4],"4th quantile",URm)
    cpm5.c2=ggHmapAPA(mats.q.c2[5],"5th quantile",URm)
    cpm6.c2=ggHmapAPA(mats.q.c2[6],"6th quantile",URm)
    cpm7.c2=ggHmapAPA(mats.q.c2[7],"7th quantile",URm)
    cpm8.c2=ggHmapAPA(mats.q.c2[8],"8th quantile",URm)
    cpm9.c2=ggHmapAPA(mats.q.c2[9],"9th quantile",URm)
    cpm10.c2=ggHmapAPA(mats.q.c2[10],"10th quantile",URm)
    
    # Measures plots
    measures.files=list.files(pattern="measures",recursive = T)
    measures.files.c1=measures.files[grep(C1,measures.files)]
    measures.files.c2=measures.files[grep(C2,measures.files)]
    
    measures.dist.c1=measures.files.c1[grep("qdist",measures.files.c1)][c(1,3:10,2)]
    measures.dist.c2=measures.files.c2[grep("qdist",measures.files.c2)][c(1,3:10,2)]
    
    measures.q=measures.files[!measures.files %in% c(measures.dist.c1,measures.dist.c2)] 
    measures.q.c1=measures.q[grep(C1,measures.q)][c(1,3:10,2)]
    measures.q.c2=measures.q[grep(C2,measures.q)][c(1,3:10,2)]
    
    md1.c1=prepMetrics(measures.dist.c1,1)
    md2.c1=prepMetrics(measures.dist.c1,2)
    md3.c1=prepMetrics(measures.dist.c1,3)
    md4.c1=prepMetrics(measures.dist.c1,4)
    md5.c1=prepMetrics(measures.dist.c1,5)
    md6.c1=prepMetrics(measures.dist.c1,6)
    md7.c1=prepMetrics(measures.dist.c1,7)
    md8.c1=prepMetrics(measures.dist.c1,8)
    md9.c1=prepMetrics(measures.dist.c1,9)
    md10.c1=prepMetrics(measures.dist.c1,10)
    md.all.c1=rbind(md1.c1,md2.c1,md3.c1,md4.c1,md5.c1,md6.c1,md7.c1,md8.c1,md9.c1,md10.c1)
    md.all.c1[,1:5]=log2(md.all.c1[,1:5])
    md.all.c1$group=C1
    
    md1.c2=prepMetrics(measures.dist.c2,1)
    md2.c2=prepMetrics(measures.dist.c2,2)
    md3.c2=prepMetrics(measures.dist.c2,3)
    md4.c2=prepMetrics(measures.dist.c2,4)
    md5.c2=prepMetrics(measures.dist.c2,5)
    md6.c2=prepMetrics(measures.dist.c2,6)
    md7.c2=prepMetrics(measures.dist.c2,7)
    md8.c2=prepMetrics(measures.dist.c2,8)
    md9.c2=prepMetrics(measures.dist.c2,9)
    md10.c2=prepMetrics(measures.dist.c2,10)
    md.all.c2=rbind(md1.c2,md2.c2,md3.c2,md4.c2,md5.c2,md6.c2,md7.c2,md8.c2,md9.c2,md10.c2)
    md.all.c2[,1:5]=log2(md.all.c2[,1:5])
    md.all.c2$group=C2
    
    md.all.merge=rbind(md.all.c1,md.all.c2)
    
    q1.c1=prepMetrics(measures.q.c1,1)
    q2.c1=prepMetrics(measures.q.c1,2)
    q3.c1=prepMetrics(measures.q.c1,3)
    q4.c1=prepMetrics(measures.q.c1,4)
    q5.c1=prepMetrics(measures.q.c1,5)
    q6.c1=prepMetrics(measures.q.c1,6)
    q7.c1=prepMetrics(measures.q.c1,7)
    q8.c1=prepMetrics(measures.q.c1,8)
    q9.c1=prepMetrics(measures.q.c1,9)
    q10.c1=prepMetrics(measures.q.c1,10)
    q.all.c1=rbind(q1.c1,q2.c1,q3.c1,q4.c1,q5.c1,q6.c1,q7.c1,q8.c1,q9.c1,q10.c1)
    q.all.c1[,1:5]=log2(q.all.c1[,1:5])
    q.all.c1$group=C1
    
    q1.c2=prepMetrics(measures.q.c2,1)
    q2.c2=prepMetrics(measures.q.c2,2)
    q3.c2=prepMetrics(measures.q.c2,3)
    q4.c2=prepMetrics(measures.q.c2,4)
    q5.c2=prepMetrics(measures.q.c2,5)
    q6.c2=prepMetrics(measures.q.c2,6)
    q7.c2=prepMetrics(measures.q.c2,7)
    q8.c2=prepMetrics(measures.q.c2,8)
    q9.c2=prepMetrics(measures.q.c2,9)
    q10.c2=prepMetrics(measures.q.c2,10)
    q.all.c2=rbind(q1.c2,q2.c2,q3.c2,q4.c2,q5.c2,q6.c2,q7.c2,q8.c2,q9.c2,q10.c2)
    q.all.c2[,1:5]=log2(q.all.c2[,1:5])
    q.all.c2$group=C2
    q.all.merge=rbind(q.all.c1,q.all.c2)
    
    p1.merge=plotMetrics(q.all.merge,!!ensym(score),xlab = "",qcpmtab)
    p2.merge=plotMetrics(md.all.merge,!!ensym(score),xlab = "",qdisttab)
    
    # report
    pdf(paste0(outdir,"/hmap_quantiles_",method,"_",score,".pdf"), width=25 ,height=6) 
    grid.arrange(arrangeGrob(cpm1.c1,cpm2.c1,cpm3.c1,cpm4.c1,cpm5.c1,cpm6.c1,cpm7.c1,cpm8.c1,cpm9.c1,cpm10.c1,nrow=1,right=C1),
                 arrangeGrob(cpm1.c2,cpm2.c2,cpm3.c2,cpm4.c2,cpm5.c2,cpm6.c2,cpm7.c2,cpm8.c2,cpm9.c2,cpm10.c2,nrow=1,right=C2),
                 p1.merge,nrow=3,top="APA analysis: by ContactCounts (cpm) quantiles")
    grid.arrange(arrangeGrob(qdist1.c1,qdist2.c1,qdist3.c1,qdist4.c1,qdist5.c1,qdist6.c1,qdist7.c1,qdist8.c1,qdist9.c1,qdist10.c1,nrow=1,right=C1),
                 arrangeGrob(qdist1.c1,qdist2.c1,qdist3.c1,qdist4.c1,qdist5.c1,qdist6.c1,qdist7.c1,qdist8.c1,qdist9.c1,qdist10.c1,nrow=1,right=C2),
                 p2.merge,nrow=3,top="APA analysis: by Loop-Distance quantiles")
    dev.off()
  }
}
