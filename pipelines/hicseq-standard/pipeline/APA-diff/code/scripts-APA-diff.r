argv = commandArgs(trailingOnly = TRUE)
inpdir = argv[1L]
APA_res = as.numeric(argv[2L])
C1 = argv[3L]
C2 = argv[4L]
URm = as.numeric(argv[5L])
outdir="plots-APA"

# inpdir="/Users/javrodher/Work/RStudio-PRJs/leukemia-cell-line-DR/data/APA-diff/APA/APA/diff/"
# outdir="/Users/javrodher/Work/RStudio-PRJs/leukemia-cell-line-DR/data/APA-diff/APA/"
# URm=2.5
# C1="CUTLL1_DMSO_A"
# C2="CUTLL1_THZ1"
# APA_res=10000

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
  print(ggplot(m,aes(x=n,y=!!ensym(metric),color=s))+
          geom_point(show.legend = F)+
          xlab(xlab)+
          ylab(ylab)+
          scale_x_continuous(breaks=c(1:nrow(m)),labels=m$s)+
          geom_hline(yintercept = log2(1),linetype=3,color="red")+
          theme(axis.text.x=element_text(size=6,angle=0,vjust = 0.5,hjust=0.5),
        plot.margin = unit(c(1, 1, 1, 1), "cm")))
}

##### RUN #####
all.files=list.files(pattern=".txt",recursive = T,full.names = F,include.dirs = F)
methods=c("APA","rankAPA","centerNormedAPA","normedAPA")
scores=c("P2M","P2UL","P2UR","P2LL","P2LR","ZscoreLL")
dirs=list.dirs(recursive = F,full.names = F)
bedpes=dirs[grep(paste0("_",C1),dirs)]
bedpes=gsub(pattern = paste0("_",C1),replacement = "",bedpes)

#method="APA"
#score="P2LL"

for (method in methods){
  print(method)
  for (score in scores){
    print(score)
    mats=all.files
    mth=paste0("^",method,".txt")
    for (bedpe in bedpes){
      
    #APA plots
    mats=list.files(pattern=mth,recursive = T)
    mats=mats[grep(bedpe,mats)]
    mat.c1=mats[grep(paste0("_",C1,"/"),mats)]
    mat.c2=mats[grep(paste0("_",C2,"/"),mats)]
    hmap.c1=ggHmapAPA(mat.c1,paste0(C1),URm = URm)
    hmap.c2=ggHmapAPA(mat.c2,paste0(C2),URm = URm)
    
    # Measures plots
    measures.files=list.files(pattern="measures",recursive = T)
    measures.files=measures.files[grep(bedpe,measures.files)]
    
    measures.file.c1=measures.files[grep(C1,measures.files)]
    measures.file.c2=measures.files[grep(C2,measures.files)]
    common.up.c1=prepMetrics(measures.file.c1,"",C1)
    common.up.c2=prepMetrics(measures.file.c2,"",C2)
   
    m.all=rbind(common.up.c1,common.up.c2)
    m.all[,1:5]=log2(m.all[,1:5])

    p1.merge=plotMetrics(m.all,!!ensym(score),xlab=bedpe,ylab = paste0("log2 (",score,")"))
    
    # report
    pdf(paste0(outdir,paste0("/hmap_diff_",bedpe,"_",method,"_",score,".pdf")), width=6 ,height=5) 
    grid.arrange(arrangeGrob(hmap.c1,hmap.c2,nrow=2,ncol=1),
                 p1.merge,nrow=1,ncol=2,top="APA analysis")
    dev.off()
    }
  }
}
