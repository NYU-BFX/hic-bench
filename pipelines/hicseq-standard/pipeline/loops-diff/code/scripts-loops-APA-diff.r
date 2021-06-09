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
read_mat = function(path) {
    mat = scan(path, "", sep = "\n", quiet = TRUE)
    mat = sapply(mat, gsub, pattern = "\\[|\\]", replacement = "")  # delete [ and ]
    mat = read.csv(text = paste(mat, collapse = "\n"), header = FALSE)
    mat = as.matrix(mat)
    colnames(mat) = NULL
    mat
}

# plot heatmap
ggHmapAPA = function(mat.path, title, nloops, URm) {
    if (file.exists(mat.path)) {
        mat = read_mat(mat.path)
        mat=rotate(mat)
        UR=mean(mat[1:6,(ncol(mat)-6):ncol(mat)])*URm
        mat = pmin(mat, UR)
        df.hmap = reshape::melt(mat)
        colnames(df.hmap) = c("row", "col", "values")
    } else {
        df.hmap=data.frame(row=1,col=1,values=0,stringsAsFactors = F)
        UR=0
    }

    axis_label = sprintf("bins (%d kb)", APA_res/1000)
    p = ggplot(df.hmap, aes(row, col, fill= values)) +
        geom_tile(show.legend = F,color=NA,size=0) +
        scale_fill_gradient(low="white", high="red", limits=c(0,UR))+
        labs(x = axis_label, y = axis_label,
             title = title, subtitle = sprintf("n=(%d)", nloops)) +
        theme(plot.title = element_text(hjust = 0.5, size=8),
              plot.subtitle = element_text(hjust = 0.5, size=8),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.line = element_blank())
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
all.files=list.files(pattern=".txt",recursive = T,full.names = F,include.dirs = F)
methods = c("APA", "rankAPA", "centerNormedAPA", "normedAPA")
scores = c("P2M", "P2UL", "P2UR", "P2LL", "P2LR", "ZscoreLL")
loopCount=read.csv(loopCount.path, sep="\t", header=F, col.names=c("count", "file"))

count.increased=loopCount$count[grep("increased",loopCount$file)]
count.decreased=loopCount$count[grep("decreased",loopCount$file)]
count.stable=loopCount$count[grep("stable",loopCount$file)]
count.specific.c1=loopCount$count[grep(paste0(C1,"_specific_loops"),loopCount$file)]
count.specific.c2=loopCount$count[grep(paste0(C2,"_specific_loops"),loopCount$file)]

measure_suffix = sprintf("%d/gw/measures.txt", APA_res)

for (method in methods){
  print(method)
  method_suffix = sprintf("%d/gw/%s.txt", APA_res, method)

  for (score in scores){
    print(score)
    mats=all.files
    mth=paste0("^",method,".txt")

    #APA plots
    mats=list.files(pattern=mth,recursive = T)

    #common
    mat.common.increased.c1=sprintf("common.loops_increased_%s/%s", C1, method_suffix)
    mat.common.decreased.c1=sprintf("common.loops_decreased_%s/%s", C1, method_suffix)
    mat.common.stable.c1=sprintf("common.loops_stable_%s/%s", C1, method_suffix)
    mat.common.increased.c2=sprintf("common.loops_increased_%s/%s", C2, method_suffix)
    mat.common.decreased.c2=sprintf("common.loops_decreased_%s/%s", C2, method_suffix)
    mat.common.stable.c2=sprintf("common.loops_stable_%s/%s", C2, method_suffix)

    common1=ggHmapAPA(mat.common.decreased.c1,"common loops - decreased",count.decreased,URm = URm)
    common2=ggHmapAPA(mat.common.decreased.c2,"common loops - decreased",count.decreased,URm = URm)
    common3=ggHmapAPA(mat.common.increased.c1,"common loops - increased",count.increased,URm = URm)
    common4=ggHmapAPA(mat.common.increased.c2,"common loops - increased",count.increased,URm = URm)
    common5=ggHmapAPA(mat.common.stable.c1,"common loops - stable",count.stable,URm = URm)
    common6=ggHmapAPA(mat.common.stable.c2,"common loops - stable",count.stable,URm = URm)

    #specific
    c1.mats.specific.c1=sprintf("%s_specific_loops_%s/%s", C1, C1, method_suffix)
    c1.mats.specific.c2=sprintf("%s_specific_loops_%s/%s", C1, C2, method_suffix)
    c2.mats.specific.c1=sprintf("%s_specific_loops_%s/%s", C2, C1, method_suffix)
    c2.mats.specific.c2=sprintf("%s_specific_loops_%s/%s", C2, C2, method_suffix)

    specific1=ggHmapAPA(c1.mats.specific.c1,sprintf("%s - specific", C1),count.specific.c1,URm = URm)
    specific2=ggHmapAPA(c1.mats.specific.c2,sprintf("%s - specific", C1),count.specific.c1,URm = URm)
    specific3=ggHmapAPA(c2.mats.specific.c1,sprintf("%s - specific", C2),count.specific.c2,URm = URm)
    specific4=ggHmapAPA(c2.mats.specific.c2,sprintf("%s - specific", C2),count.specific.c2,URm = URm)

    # Measures plots
    measures.files=list.files(pattern="measures",recursive = T)
    measures.files=measures.files[grep("increased|decreased|stable|specific",measures.files)]

    #common
    mat.common.increased.c1=sprintf("common.loops_increased_%s/%s", C1, measure_suffix)
    mat.common.decreased.c1=sprintf("common.loops_decreased_%s/%s", C1, measure_suffix)
    mat.common.stable.c1=sprintf("common.loops_stable_%s/%s", C1, measure_suffix)
    mat.common.increased.c2=sprintf("common.loops_increased_%s/%s", C2, measure_suffix)
    mat.common.decreased.c2=sprintf("common.loops_decreased_%s/%s", C2, measure_suffix)
    mat.common.stable.c2=sprintf("common.loops_stable_%s/%s", C2, measure_suffix)

    common.up.c1=prepMetrics(mat.common.increased.c1,"common increased",C1)
    common.up.c2=prepMetrics(mat.common.increased.c2,"common increased",C2)
    common.down.c1=prepMetrics(mat.common.decreased.c1,"common decreased",C1)
    common.down.c2=prepMetrics(mat.common.decreased.c2,"common decreased",C2)
    common.stable.c1=prepMetrics(mat.common.stable.c1,"common stable",C1)
    common.stable.c2=prepMetrics(mat.common.stable.c2,"common stable",C2)
    common.merged=rbind(common.down.c1,common.down.c2,common.up.c1,common.up.c2,common.stable.c1,common.stable.c2)

    #specific
    c1.measures.files.specific.c1=sprintf("%s_specific_loops_%s/%s", C1, C1, measure_suffix)
    c1.measures.files.specific.c2=sprintf("%s_specific_loops_%s/%s", C1, C2, measure_suffix)
    c2.measures.files.specific.c1=sprintf("%s_specific_loops_%s/%s", C2, C1, measure_suffix)
    c2.measures.files.specific.c2=sprintf("%s_specific_loops_%s/%s", C2, C2, measure_suffix)

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
