argv = commandArgs(trailingOnly = TRUE)
inpdir = argv[1L]
APA_res = as.numeric(argv[2L])
C1 = argv[3L]
C2 = argv[4L]
URm = as.numeric(argv[5L])
loopCount.path = argv[6L]
outdir="in-house-analysis"

library(ggplot2)
library(gridExtra)
library(pheatmap)
library(RColorBrewer)
library(grid)
library(ggplotify)

setwd(inpdir)
options(scipen=10000)
dir.create(outdir)
pdf(NULL)

#### Functions ####

# rotate matrix
rotate=function(x) {t(apply(x, 2, rev))}

# plot heatmap
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

# prepare scores table
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

# plot scores by loop-subset
plotMetrics=function(m,metric,xlab,ylab){
  m$n=1:nrow(m)
  print(ggplot(m,aes(x=n,y=!!ensym(metric),color=q))+
          geom_point(show.legend = F)+
          xlab(xlab)+
          ylab(ylab)+
          scale_x_continuous(breaks=c(1:nrow(m)),labels=m$q,sec.axis = dup_axis(labels=m$s))+
          geom_hline(yintercept = log2(1),linetype=3,color="red")+
          theme(axis.text.x=element_text(size=6,angle=0,vjust = 0.5,hjust=0.5)))
}

##### RUN #####

# set labels (to avoid ggplot error when first character of label is numeric -> add "S" as first char)
test.1st.char=as.numeric(substr(x = C1,start = 1,stop = 1)) 
test.1st.char=test.1st.char[!is.na(test.1st.char)]
if(length(test.1st.char)>0){C1=paste0("S",C1)}
test.1st.char=as.numeric(substr(x = C2,start = 1,stop = 1)) 
test.1st.char=test.1st.char[!is.na(test.1st.char)]
if(length(test.1st.char)>0){C2=paste0("S",C2)}

# load data
all.files=list.files(pattern=".txt",recursive = T,full.names = F,include.dirs = F)
methods=c("APA","rankAPA","centerNormedAPA","normedAPA")
scores=c("P2M","P2UL","P2UR","P2LL","P2LR","ZscoreLL")
loopCount=read.csv(loopCount.path,sep="\t",header=F)
names(loopCount)=c("count","file")
count.increased=loopCount$count[grep("increased",loopCount$file)]
count.decreased=loopCount$count[grep("decreased",loopCount$file)]
count.stable=loopCount$count[grep("stable",loopCount$file)]
count.specific.c1=loopCount$count[grep(paste0(C1,"_specific_loops"),loopCount$file)]
count.specific.c2=loopCount$count[grep(paste0(C2,"_specific_loops"),loopCount$file)]

for (method in methods){
  print(method)
  for (score in scores){
    print(score)
    mats=all.files
    mth=paste0("^",method,".txt")
    
    #APA plots
    mats=list.files(pattern=mth,recursive = T)
    
    #common
    mats.common=mats[grep("common",mats)]
    mats.common=mats.common[grep("increased|decreased|stable",mats.common)]
    mats.common.c1=mats.common[grep(C1,mats.common)]
    mats.common.c2=mats.common[grep(C2,mats.common)]
    
    mat.common.increased.c1=mats.common.c1[grep("increased",mats.common.c1)]
    mat.common.decreased.c1=mats.common.c1[grep("decreased",mats.common.c1)]
    mat.common.stable.c1=mats.common.c1[grep("stable",mats.common.c1)]
    mat.common.increased.c2=mats.common.c2[grep("increased",mats.common.c2)]
    mat.common.decreased.c2=mats.common.c2[grep("decreased",mats.common.c2)]
    mat.common.stable.c2=mats.common.c2[grep("stable",mats.common.c2)]
    
    common1=ggHmapAPA(mat.common.decreased.c1,paste0("common loops - decreased \n (n=",count.decreased,")"),URm = URm)
    common2=ggHmapAPA(mat.common.decreased.c2,paste0("common loops - decreased \n (n=",count.decreased,")"),URm = URm)
    common3=ggHmapAPA(mat.common.increased.c1,paste0("common loops - increased \n (n=",count.increased,")"),URm = URm)
    common4=ggHmapAPA(mat.common.increased.c2,paste0("common loops - increased \n (n=",count.increased,")"),URm = URm)
    common5=ggHmapAPA(mat.common.stable.c1,paste0("common loops - stable \n (n=",count.stable,")"),URm = URm)
    common6=ggHmapAPA(mat.common.stable.c2,paste0("common loops - stable \n (n=",count.stable,")"),URm = URm)
    
    #specific
    mats.specific=mats[!mats %in% mats.common]
    c1.mats.specific=mats.specific[grep(paste0(C1,"_specific_loops"),mats.specific)]
    c2.mats.specific=mats.specific[grep(paste0(C2,"_specific_loops"),mats.specific)]
    
    c1.mats.specific.c1=c1.mats.specific[grep(paste0("_",C1,"/"),c1.mats.specific,fixed = T)]
    c1.mats.specific.c2=c1.mats.specific[grep(paste0("_",C2,"/"),c1.mats.specific,fixed = T)]
    c2.mats.specific.c1=c2.mats.specific[grep(paste0("_",C1,"/"),c2.mats.specific,fixed = T)]
    c2.mats.specific.c2=c2.mats.specific[grep(paste0("_",C2,"/"),c2.mats.specific,fixed = T)]
    
    specific1=ggHmapAPA(c1.mats.specific.c1,paste0(C1,"\nspecific\n (n=",count.specific.c1,")"),URm = URm)
    specific2=ggHmapAPA(c1.mats.specific.c2,paste0(C1,"\nspecific\n (n=",count.specific.c1,")"),URm = URm)
    specific3=ggHmapAPA(c2.mats.specific.c1,paste0(C2,"\nspecific\n (n=",count.specific.c2,")"),URm = URm)
    specific4=ggHmapAPA(c2.mats.specific.c2,paste0(C2,"\nspecific\n (n=",count.specific.c2,")"),URm = URm)
    
    # Measures plots
    measures.files=list.files(pattern="measures",recursive = T)
    measures.files=measures.files[grep("increased|decreased|stable|specific",measures.files)]
    
    #common
    measures.files.common=measures.files[grep("common",measures.files)]
    measures.files.common.c1=measures.files.common[grep(C1,measures.files.common)]
    measures.files.common.c2=measures.files.common[grep(C2,measures.files.common)]
    
    mat.common.increased.c1=measures.files.common.c1[grep("increased",measures.files.common.c1)]
    mat.common.decreased.c1=measures.files.common.c1[grep("decreased",measures.files.common.c1)]
    mat.common.stable.c1=measures.files.common.c1[grep("stable",measures.files.common.c1)]
    mat.common.increased.c2=measures.files.common.c2[grep("increased",measures.files.common.c2)]
    mat.common.decreased.c2=measures.files.common.c2[grep("decreased",measures.files.common.c2)]
    mat.common.stable.c2=measures.files.common.c2[grep("stable",measures.files.common.c2)]
    
    common.up.c1=prepMetrics(mat.common.increased.c1,"common increased",C1)
    common.up.c2=prepMetrics(mat.common.increased.c2,"common increased",C2)
    common.down.c1=prepMetrics(mat.common.decreased.c1,"common decreased",C1)
    common.down.c2=prepMetrics(mat.common.decreased.c2,"common decreased",C2)
    common.stable.c1=prepMetrics(mat.common.stable.c1,"common stable",C1)
    common.stable.c2=prepMetrics(mat.common.stable.c2,"common stable",C2)
    common.merged=rbind(common.down.c1,common.down.c2,common.up.c1,common.up.c2,common.stable.c1,common.stable.c2)
    
    #specific
    measures.files.specific=measures.files[!measures.files %in% measures.files.common]
    c1.measures.files.specific=measures.files.specific[grep(paste0(C1,"_specific_loops"),measures.files.specific)]
    c2.measures.files.specific=measures.files.specific[grep(paste0(C2,"_specific_loops"),measures.files.specific)]
    
    c1.measures.files.specific.c1=c1.measures.files.specific[grep(paste0("_",C1,"/"),c1.measures.files.specific,fixed = T)]
    c1.measures.files.specific.c2=c1.measures.files.specific[grep(paste0("_",C2,"/"),c1.measures.files.specific,fixed = T)]
    c2.measures.files.specific.c1=c2.measures.files.specific[grep(paste0("_",C1,"/"),c2.measures.files.specific,fixed = T)]
    c2.measures.files.specific.c2=c2.measures.files.specific[grep(paste0("_",C2,"/"),c2.measures.files.specific,fixed = T)]
    
    md1=prepMetrics(c1.measures.files.specific.c1,paste0(C1,"\nspecific"),C1)
    md2=prepMetrics(c1.measures.files.specific.c2,paste0(C1,"\nspecific"),C2)
    md3=prepMetrics(c2.measures.files.specific.c1,paste0(C2,"\nspecific"),C1)
    md4=prepMetrics(c2.measures.files.specific.c2,paste0(C2,"\nspecific"),C2)
    md.all=rbind(md1,md2,md3,md4)
    
    common.specific=rbind(common.merged,md.all)
    p1=plotMetrics(common.specific,!!ensym(score),xlab = "",ylab = paste0("log2 (",score,")"))
    
    # report
    pdf(paste0(outdir,paste0("/hmap_diff_",method,"_",score,".pdf")), width=10 ,height=6) 
    grid.arrange(arrangeGrob(common1,common3,common5,specific1,specific3,nrow=1,left=C1,right=C1),
                 arrangeGrob(common2,common4,common5,specific2,specific4,nrow=1,left=C2,right=C2),
                 p1,nrow=3,top="APA analysis on FitHiC loop-subsets")
    dev.off()
  }
}
