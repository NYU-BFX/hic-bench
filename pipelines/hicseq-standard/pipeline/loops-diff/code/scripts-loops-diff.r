argv = commandArgs(trailingOnly = T)
loops.C1 = argv[1L]
loops.C2 = argv[2L]
C1.lab = argv[3L]
C2.lab = argv[4L]
outdir = argv[5L]
binsize = as.numeric(argv[6L])
tss.path = argv[7L]
common.lfc = as.numeric(argv[8L])
qcut1 = as.numeric(argv[9L])
qcut2 = as.numeric(argv[10L])
min.dist = as.numeric(argv[11L])
max.dist = as.numeric(argv[12L])
annot_promoters = F # to do

library(ggplot2)
library(plyr)
library(ggsignif)
library(gridExtra)
library(GenomicRanges)

### Functions ##############################################################
pdf(NULL)
options(scipen=10000)
# plot loop contactCounts by distance
contactsByDist=function(df.gg,group,is.common){
  x=group
  if(is.common){
    color.df=keys.df.common
    title="Common loops"
  } else { if(!!ensym(x)=="group"){
    color.df=keys.df
    title="All loops (by sample/group)"
  } else {
    color.df=keys.df2
    title="Common & group-specific loops"
  } 
  }
  
  print(ggplot(df.gg, aes(x=distance, y=contactCount,group=!!ensym(x),color=!!ensym(x)))+
          scale_x_continuous(trans='log10')+
          scale_y_continuous(trans= "log10")+
          geom_smooth(size=0.6,se = T)+
          ggtitle(label=title)+
          scale_color_manual(values = color.df$color)+
          xlab("distance")+
          ylab("contactCount")+
          xlim(0,2500000)+
          theme_bw() +  theme(plot.title = element_text(hjust = 0.5,size=12),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank()) +
          theme(axis.text.x = element_text(angle = 90),legend.title = element_blank()))
}

# plot loop contactCounts by distance by chromosome
contactsByDistByChr=function(df.gg,group,is.common){
  x=group
  if(is.common){
    color.df=keys.df.common
    title="Common loops"
  } else { if(!!ensym(x)=="group"){
    color.df=keys.df
    title="All loops (by sample/group)"
  } else {
    color.df=keys.df2
    title="Common & group-specific loops"
  } 
  }
  
  print(ggplot(df.gg, aes(x=distance, y=contactCount,group=!!ensym(x),color=!!ensym(x)))+
          scale_x_continuous(trans='log10')+
          scale_y_continuous(trans= "log10")+
          geom_smooth(size=0.6)+
          ggtitle(label=title)+
          scale_color_manual(values = color.df$color)+
          facet_wrap(chrom~., scales= "fixed",nrow=3)+
          xlab("distance")+
          ylab("contactCount")+
          xlim(0,2500000)+
          theme_bw() +  theme(plot.title = element_text(hjust = 0.5,size=12),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank()) +
          theme(axis.text.x = element_text(angle = 90),legend.title = element_blank()))
}

# hexbin plot of contactCounts
hexbinContacts=function(df){
  correl=round(cor(df.common$C1,df.common$C2),digits = 3)
  print(ggplot(data = df.common, aes_string(x="C1",y="C2"))+ 
          stat_binhex(bins=70) +
          ggtitle("Common loops")+
          # scale_x_continuous(trans='log10')+
          #scale_y_continuous(trans= "log10")+
          xlab(paste0("contactCount [",C1.lab,"]"))+
          ylab(paste0("contactCount [",C2.lab,"]"))+
          scale_fill_gradientn(colours=c("blue","orange","red"),trans="log10","  log10 count\n(contactCounts)\n")+
          geom_abline(slope = 1,intercept = 0,linetype=2)+
          annotate(geom = 'text', label = paste0('r = ',correl), x = -Inf, y = Inf, hjust = -0.2, vjust = 1.5,size=3)+
          theme(plot.title = element_text(color = "black",hjust = 0.5),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                axis.line = element_line(colour = "black")))
}

# boxplot loop distance
distanceBox=function(df.common, max.dist){
  df.gg.distance=data.frame(distance=c(df.common$distance,C1.specific$distance,C2.specific$distance),group=c(rep("common",nrow(df.common)),rep(C1.lab,nrow(C1.specific)),rep(C2.lab,nrow(C2.specific))),stringsAsFactors = F)
  df.gg.distance=df.gg.distance[df.gg.distance$distance < max.dist,] #filter loops < 10Mb
  df.gg.distance$group = as.factor(df.gg.distance$group)
  
  p_meds = ddply(df.gg.distance, .(group), summarise, med = median(distance,na.rm=T))
  p_meds$med = round(p_meds$med,digits = 3)
  p_meds$med=p_meds$med/1000
  if(nrow(p_meds)>1){tests=split(t(combn(levels(df.gg.distance$group), 2)), seq(nrow(t(combn(levels(df.gg.distance$group), 2)))))} else{tests=NA}
  min.y=min(df.gg.distance$distance,na.rm=T)
  max.y=max(df.gg.distance$distance,na.rm=T)
  title="Loop distance in common and group-specific loops"
  subtitle=paste0("Median loop distance (kb): Common=",p_meds$med[1],"; C2-specific=",p_meds$med[2],"; C1-specific=",p_meds$med[3])
  
  #set colors
  df.gg.distance$group=as.factor(df.gg.distance$group)
  df.gg.distance$group=factor(df.gg.distance$group,levels(df.gg.distance$group))
  ins.stat=unique(as.character(df.gg.distance$group))
  df.gg.distance$group=factor(df.gg.distance$group,c("common",ins.stat[ins.stat != "common"]))
  
  keys.df.distance=data.frame(keys=levels(df.gg.distance$group),color=NA)
  keys.df.distance$color[keys.df.distance$keys =="common"]="darkgrey"
  keys.df.distance$color[keys.df.distance$keys == C1.lab]="blue"
  keys.df.distance$color[keys.df.distance$keys == C2.lab]="red"
  
  print(ggplot(df.gg.distance, aes(x=group, y=distance, group=group, fill=group))+
          xlab(label = "") + 
          ylab(label = "loop distance (bp)") + 
          ylim(c(min.y,max.y))+
          scale_fill_manual(values = keys.df.distance$color)+
          labs(fill="loop class")+
          ggtitle(label=title,subtitle = subtitle )+
          geom_boxplot(width=0.5,outlier.size = 0.01)+
          theme(plot.title = element_text(hjust = 0.5,size=12),
                plot.subtitle = element_text(color = "black",hjust = 0.5,size = 8.5),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                panel.background = element_blank(), axis.line = element_line(colour = "black"))+
          theme(axis.text.x = element_text(angle = -45))+
          geom_signif(comparisons = tests, map_signif_level=TRUE,test = 'wilcox.test',na.rm = T,y_position = max.y*1))
}

# loop count by distance
loopCountByDist=function(df){
  df=df.distance.num
  #set colors
  df$group=as.factor(df$group)
  df$group=factor(df$group,levels(df$group))
  ins.stat=unique(as.character(df$group))
  df$group=factor(df$group,c("common",ins.stat[ins.stat != "common"]))
  keys.df.distance=data.frame(keys=levels(df$group),color=NA)
  keys.df.distance$color[keys.df.distance$keys =="common"]="darkgrey"
  keys.df.distance$color[keys.df.distance$keys == C1.lab]="blue"
  keys.df.distance$color[keys.df.distance$keys == C2.lab]="red"
  
  print(ggplot(df, aes(x=distance, y=count,group=group,color=group))+
          scale_x_continuous(trans='log10')+
          scale_y_continuous(trans= "log10")+
          geom_smooth(size=0.3)+
          scale_color_manual(name="loop class",values = keys.df.distance$color)+
          xlab("distance")+
          ylab("loop count")+
          xlim(0,5000000)+
          labs()+
          theme_bw() +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank()) +
          theme(axis.text.x = element_text(angle = 90)))
}

# lfc contacts by rank
lfcRankedContacts = function(df.common){
  ggplot(data=df.common,aes(x=lfc.contacts,y=rank))+
    geom_smooth()+
    #facet_wrap(chr1~., scales= "fixed",ncol=4)
    geom_vline(xintercept = 0,linetype="dotted")+
    geom_vline(xintercept = 1,linetype="dotted",color="red")+
    geom_vline(xintercept = -1,linetype="dotted",color="red")+
    ylab(paste0("common loops ranked (",C1.lab,")"))+
    xlab(paste0("log2FC contactCounts (",C2.lab,"/",C1.lab,")"))
}

# scatter plot contactCount
scatterContactCount = function(df.common){
  print(ggplot(df.common, aes(x=C1, y=C2,color=loop.strength))+ 
          geom_point(size=2)+
          #scale_x_continuous(trans='log10')+
          #scale_y_continuous(trans= "log10")+
          ggtitle(label="Common loops",subtitle=paste0("log2FC(contactCounts) threshold = ",common.lfc))+
          scale_color_manual(values=keys.df.common2$color)+
          geom_abline(slope = 1,intercept = 0,linetype=3,color="black")+
          xlab(paste0("contactCount [",C1.lab,"]"))+
          ylab(paste0("contactCount [",C2.lab,"]"))+
          theme(plot.title = element_text(color = "black",hjust = 0.5),plot.subtitle = element_text(color = "black",hjust = 0.5,size = 8.5),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")))
}

# loop count barplot (common, C1.specific, C2.specific)
countBarLoopClass = function(df.gg.bar){
  print(ggplot(data=df.gg.bar,aes(x=group2,fill=group2))+
          geom_bar(colour="black")+
          scale_fill_manual(values = df.count$color)+
          xlab("")+
          ggtitle("")+
          ylab("loop count")+
          labs(fill="loop class")+
          theme(plot.title = element_text(color = "black",hjust = 0.5),axis.text.x = element_text(angle = -45))+
          geom_text(data=df.count,aes(x=group2,y=count,label=count),size=3,vjust=1.4,color="white"))
}

# loop count barplot (common, C1.specific, C2.specific)
countBarStacked = function(df.gg2){
  print(ggplot(data=df.gg2,aes(x=group,fill=group2))+
          geom_bar(colour="black")+
          ylab("loop count")+
          xlab("")+
          labs(fill="loop class")+
          theme(axis.text.x = element_text(angle = -45))+
          scale_fill_manual(values = keys.df3$color))
}

# loop % barplot (common, C1.specific, C2.specific)
fractionBarStacked = function(df.perc){
  print(ggplot(data=df.perc,aes(x=group,fill=group2,y=perc))+
          geom_bar(position="fill",stat="identity",colour="black")+
          ylab("loop fraction")+
          xlab("")+
          labs(fill="loop class")+
          scale_fill_manual(values = c("darkgrey","purple"))+
          coord_flip())
}

# loop % barplot (common, C1.specific, C2.specific)
fractionBarStackedLoopClass = function(df.count){
  print(ggplot(data=df.count,aes(x=all,fill=group2,y=perc))+
          geom_bar(position="fill",stat="identity",colour="black",width=0.45)+
          scale_fill_manual(values = df.count$color)+
          xlab("")+
          labs(fill="loop class")+
          ylab("loop fraction")+
          coord_flip())
}

# annotate loops-promoter 
annotGenes=function(x){
  loopIndex=x
  print(x)
  TssIndex.1 = anchor1.tss.index$subjectHits[anchor1.tss.index$queryHits==loopIndex]
  TssIndex.2 = anchor2.tss.index$subjectHits[anchor2.tss.index$queryHits==loopIndex]
  current.annot.1 = df.gg$anchor1.annot[loopIndex]
  current.annot.2 = df.gg$anchor2.annot[loopIndex]
  
  #anchor1
  if(df.gg$anchor1.promoter[loopIndex]==T){
    if(!is.na(current.annot.1)){ 
      A1.annot=paste(current.annot.1,tss$geneID[TssIndex.1],collapse = ",")
    } else{
      A1.annot=paste(tss$geneID[TssIndex.1],collapse = ",")
    }
  } else {A1.annot=current.annot.1}
  #anchor2
  if(df.gg$anchor2.promoter[loopIndex]==T){
    if(!is.na(current.annot.2)){ 
      A2.annot=paste(current.annot.2,tss$geneID[TssIndex.2],collapse = ",")
    } else{
      A2.annot=paste(tss$geneID[TssIndex.2],collapse = ",")
    }
  } else {A2.annot=current.annot.2}
  return(list(a1=A1.annot,a2=A2.annot))
}

# promoter count
promoterCount = function(df.gg.bar){
  ggplot(data=df.gg.bar,aes(x=group2,fill=any.promoter))+
    geom_bar(colour="black",width=0.45)+
    labs(fill="")+
    ylab("loop count")+
    theme(axis.text.x = element_text(angle = -45))+
    xlab("")
}

# promoter fraction
promoterFraction = function(df.promoter.perc){
  ggplot(data=df.promoter.perc,aes(x=group2,y=frac,fill=group))+
    geom_bar(position="fill",stat="identity",colour="black",width=0.45)+
    labs(fill="")+
    ylab("loop fraction")+
    theme(axis.text.x = element_text(angle = -45))+
    xlab("")
}

# promoter count
promoterCountCommon = function(df.common){
  df.common.nostable=df.common[df.common$loop.strength != "stable",]
  ggplot(data=df.common.nostable,aes(x=loop.strength,fill=any.promoter))+
    geom_bar(colour="black",width=0.45)+
    labs(fill="")+
    ylab("loop count")+
    labs(fill="common loops")+
    theme(axis.text.x = element_text(angle = -45))+
    xlab("")
}

# promoter fraction common
promoterFractionCommon = function(df.promoter.perc.common){
  ggplot(data=df.promoter.perc.common,aes(x=group2,y=frac,fill=group))+
    geom_bar(position="fill",stat="identity",colour="black",width=0.45)+
    labs(fill="")+
    ylab("loop fraction")+
    labs(fill="common loops")+
    theme(axis.text.x = element_text(angle = -45))+
    xlab("")
}

### START ###################################################
dir.create(paste0(outdir,"/plots/"))

### Prepare data ###
# load loop data
C1.loops=read.table(loops.C1,header = T,stringsAsFactors = F)
C2.loops=read.table(loops.C2,header = T, stringsAsFactors = F)
tss = read.table(tss.path,header=F,stringsAsFactors = F)
tss = tss[,1:4]
colnames(tss)=c("chr","start","end","geneID")
tss$geneID=gsub(pattern = "tss:",replacement = "",x = tss$geneID)

# set labels (to avoid ggplot error when first character of label is numeric -> add "S" as first char)
C1.lab.original=C1.lab
C2.lab.original=C2.lab
test.1st.char=as.numeric(substr(x = C1.lab,start = 1,stop = 1)) 
test.1st.char=test.1st.char[!is.na(test.1st.char)]
if(length(test.1st.char)>0){C1.lab=paste0("S",C1.lab)}
test.1st.char=as.numeric(substr(x = C2.lab,start = 1,stop = 1)) 
test.1st.char=test.1st.char[!is.na(test.1st.char)]
if(length(test.1st.char)>0){C2.lab=paste0("S",C2.lab)}

# generate master table
C1.loops$distance=abs(C1.loops$fragmentMid2-C1.loops$fragmentMid1)
C2.loops$distance=abs(C2.loops$fragmentMid2-C2.loops$fragmentMid1)
C1.loops$group=C1.lab
C2.loops$group=C2.lab
C1.loops$loops.ID=paste0(C1.loops$chr1,":",C1.loops$fragmentMid1,":",C1.loops$fragmentMid2)
C2.loops$loops.ID=paste0(C2.loops$chr1,":",C2.loops$fragmentMid1,":",C2.loops$fragmentMid2)

df.gg=rbind(C1.loops,C2.loops)
df.gg=df.gg[df.gg$distance >= min.dist & df.gg$distance <= max.dist,] #filter by distance cutoffs
df.gg$chr=df.gg$chr1
df.gg$contactCount=as.numeric(df.gg$contactCount)
df.gg$distance=as.numeric(df.gg$distance)
df.gg$group=as.factor(df.gg$group)
df.gg$chrom=as.numeric(gsub(df.gg$chr,pattern = "chr",replacement = ""))
df.gg=df.gg[order(df.gg$chrom),]
df.gg$chrom[is.na(df.gg$chrom)]="X"
df.gg$chrom=paste0("chr",df.gg$chrom)
chrom=unique(df.gg$chrom)
df.gg$chrom=as.factor(df.gg$chrom)
df.gg$chrom=factor(df.gg$chrom,levels=chrom)
df.gg$loops.ID=paste0(df.gg$chr1,":",df.gg$fragmentMid1,":",df.gg$fragmentMid2)

# Annotate promoters
tss.gr=makeGRangesFromDataFrame(tss,keep.extra.columns = T)
df.gg$start1=df.gg$fragmentMid1
df.gg$end1=df.gg$fragmentMid1
df.gg$start2=df.gg$fragmentMid2
df.gg$end2=df.gg$fragmentMid2
df.gg$anchor1.promoter=F
df.gg$anchor2.promoter=F
df.gg$anchor1.annot=NA
df.gg$anchor2.annot=NA

anchor1.gr=makeGRangesFromDataFrame(df.gg[,c(12,13,15:16)],keep.extra.columns = T,start.field = "start1",end.field = "end1")
anchor2.gr=makeGRangesFromDataFrame(df.gg[,c(12,13,17:18)],keep.extra.columns = T,start.field = "start2",end.field = "end2")
anchor1.tss.index=as.data.frame(findOverlaps(anchor1.gr,tss.gr,maxgap = binsize/2))
anchor2.tss.index=as.data.frame(findOverlaps(anchor2.gr,tss.gr,maxgap = binsize/2))
df.gg$anchor1.promoter[anchor1.tss.index$queryHits]=T
df.gg$anchor2.promoter[anchor2.tss.index$queryHits]=T

if (annot_promoters){ ## This takes too long: Look for a different approach.
  loops.A1.p = df.gg[df.gg$anchor1.promoter==T,]
  loops.A2.p = df.gg[df.gg$anchor2.promoter==T,]
  loops.index=as.numeric(rownames(df.gg[df.gg$anchor1.promoter==T | df.gg$anchor2.promoter==T,]))
  for (i in loops.index){ df.gg[i,21:22]=unlist(annotGenes(x = i)) }
}

# classify loops in common & group-specific
df.gg$qcut1=F
df.gg$qcut1[df.gg$q.value <= qcut1]=T 
C1.loops$qcut1=F
C1.loops$qcut1[C1.loops$q.value <= qcut1]=T 
C2.loops$qcut1=F
C2.loops$qcut1[C2.loops$q.value <= qcut1]=T 

x=C1.loops$loops.ID[C1.loops$loops.ID %in% C2.loops$loops.ID]
c1=C1.loops[C1.loops$loops.ID %in% x,]
c2=C2.loops[C2.loops$loops.ID %in% x,]
c1=c1[order(c1$loops.ID),]
c2=c2[order(c2$loops.ID),]
c1$c2.qcut1=c2$qcut1
common.loops=c1$loops.ID[c1$qcut1==T | c1$c2.qcut1==T]

# Criteria: e.g. group1-specific -> L1 loops conditions: not in common.loops, in L1 set with qval < qcut1, not in L2 with qval < qcut1
y=df.gg[!(df.gg$loops.ID %in% common.loops),]
C1.specific=y[y$loops.ID %in% C1.loops$loops.ID[C1.loops$qcut1==T] & !(y$loops.ID %in% C2.loops$loops.ID[C2.loops$qcut1==T]),]
C1.specific=C1.specific[C1.specific$qcut1==T,]
C2.specific=y[y$loops.ID %in% C2.loops$loops.ID[C2.loops$qcut1==T] & !(y$loops.ID %in% C1.loops$loops.ID[C1.loops$qcut1==T]),]
C2.specific=C2.specific[C2.specific$qcut1==T,]

df.gg=df.gg[df.gg$loops.ID %in% c(C1.specific$loops.ID,C2.specific$loops.ID,common.loops),]
df.gg$group2=as.character(df.gg$group)
df.gg$group2[df.gg$loops.ID %in% common.loops]="common"
df.gg$group2[df.gg$group2 == "common" & df.gg$group == C1.lab]=paste0("common-",C1.lab)
df.gg$group2[df.gg$group2 == "common" & df.gg$group == C2.lab]=paste0("common-",C2.lab)

master.tab=df.gg[df.gg$group2 == "common" | df.gg$qcut1 == T,]
master.tab$group2=as.factor(master.tab$group2)
master.tab=master.tab[,-c(23)]
df.gg=df.gg[,-c(23)]

write.table(C1.specific[c(1:5,7,10,19:20)],paste0(outdir,"/",C1.lab.original,"_specific_loops.tsv"),row.names = F,col.names = T,quote = F,sep="\t")
write.table(C2.specific[c(1:5,7,10,19:20)],paste0(outdir,"/",C2.lab.original,"_specific_loops.tsv"),row.names = F,col.names = T,quote = F,sep="\t")
write.table(master.tab[,c(1:5,7,10,11,19:20,23)],paste0(outdir,"/master_table.tsv"),row.names = F,col.names = T,quote = F,sep="\t")

# common loops
df.gg.common=df.gg[df.gg$loops.ID %in% common.loops,]
common.C1=df.gg.common[df.gg.common$group==C1.lab,]
common.C2=df.gg.common[df.gg.common$group==C2.lab,]
common.C1$C1=common.C1$contactCount
common.C2$C2=common.C2$contactCount

df.common = cbind(common.C1,common.C2)
df.common=df.common[,c(1:4,11,24,34,48:47,43:44)]
df.common$loop.strength=NA
df.common$lfc.contacts=log2(df.common$C2/df.common$C1)
df.common$loop.strength[df.common$lfc.contacts >= common.lfc]="increased"
df.common$loop.strength[df.common$lfc.contacts <= -common.lfc]="decreased"
df.common$loop.strength[abs(df.common$lfc.contacts) < common.lfc]="stable"
df.common=df.common[order(df.common$C1,decreasing = T),]
df.common$rank=nrow(df.common):1

common.stable=df.common[df.common$loop.strength=="stable",c(1:4,7,10:11,13)]
common.up=df.common[df.common$loop.strength=="increased",c(1:4,7,10:11,13)]
common.down=df.common[df.common$loop.strength=="decreased",c(1:4,7,10:11,13)]
write.table(common.up,paste0(outdir,"/common.loops_increased.tsv"),row.names = F,col.names = T,quote = F,sep="\t")
write.table(common.down,paste0(outdir,"/common.loops_decreased.tsv"),row.names = F,col.names = T,quote = F,sep="\t")
write.table(common.stable,paste0(outdir,"/common.loops_stable.tsv"),row.names = F,col.names = T,quote = F,sep="\t")
write.table(df.common,paste0(outdir,"/common.loops.tsv"),row.names = F,col.names = T,quote = F,sep="\t")

#distance
df.distance.count=data.frame(distance=unique(df.gg$distance),count.C1=NA,count.C2=NA,count.common=NA,stringsAsFactors = F)
df.distance.count=df.distance.count[order(df.distance.count$distance,decreasing = F),]
df.distance.count=df.distance.count[df.distance.count$distance < max.dist,]

for (dist.num in 1:nrow(df.distance.count)){
  catch=df.gg[df.gg$distance == df.distance.count$distance[dist.num],]
  count.specific.C1=sum(catch$group2 == C1.lab,na.rm=T)
  count.specific.C2=sum(catch$group2 == C2.lab,na.rm=T)
  count.common=sum(catch$group2 == "common",na.rm=T)
  df.distance.count$count.C1[dist.num]=count.specific.C1
  df.distance.count$count.C2[dist.num]=count.specific.C2
  df.distance.count$count.common[dist.num]=count.common
}
df.distance.num=data.frame(distance=rep(df.distance.count$distance,3),count=c(df.distance.count$count.C1,df.distance.count$count.C2,df.distance.count$count.common),group=c(rep(C1.lab,nrow(df.distance.count)),rep(C2.lab,nrow(df.distance.count)),rep("common",nrow(df.distance.count))),stringsAsFactors = F)

# barplot data frame 
df.gg.bar=df.gg
c1.unique.count=sum(df.gg.bar$group2==C1.lab,na.rm = T)
c2.unique.count=sum(df.gg.bar$group2==C2.lab,na.rm = T)
df.gg.bar=df.gg.bar[!duplicated(df.gg.bar$loops.ID),]
common.count=length(grep("common",df.gg.bar$group2))

# get counts
df.gg.bar2=df.gg.bar
df.gg.bar2$group2[grep("common",df.gg.bar2$group2)]="common"
df.count=data.frame(group2=c("common",C1.lab,C2.lab),count=c(common.count,c1.unique.count,c2.unique.count))
df.count$group2=as.factor(df.count$group2)
df.count$group2=factor(df.count$group2,levels=c("common",C1.lab,C2.lab))
df.count$color=c("darkgrey","blue","red")
df.count=df.count[df.count$group2 %in% df.gg.bar2$group2,]
tot.count=sum(df.count$count,na.rm = T)
df.count$perc=round(df.count$count/tot.count,digits=2)
df.count$all="loop class"
df.gg.bar2$group2=factor(df.gg.bar2$group2,levels = df.count$group2)
df.gg.bar2$all="loop class"
df.count$comparison=paste0(C2.lab,"_vs_",C1.lab)
#write.table(df.count,paste0(outdir,"/metrics_count.tsv"),row.names = F,col.names = T,quote = F,sep="\t")

# get fractions
df.perc=data.frame(group2=c(C1.lab,C2.lab),unique=NA,common=NA)
x=df.count$group2[grep("common",df.count$group2,invert = T)]
for (i in x){
  df.perc$unique[df.perc$group2 == i]=round(df.count$count[df.count$group2 == i]/(df.count$count[df.count$group2 == i]+df.count$count[grep("common",df.count$group2)]),digits = 2)
  df.perc$common[df.perc$group2 == i]=1-df.perc$unique[df.perc$group2 == i]
}

df.perc=data.frame(group=rep(df.perc$group2,2),perc=c(df.perc$unique,df.perc$common))
df.perc$group2=rep(c("common","unique"),each=2)
df.perc$perc=df.perc$perc*100
df.perc$comparison=paste0(C2.lab,"_vs_",C1.lab)
#write.table(df.perc,paste0(outdir,"/metrics_fraction.tsv"),row.names = F,col.names = T,quote = F,sep="\t")

# promoter count / fraction
df.gg.bar2$any.promoter="no-promoter"
df.gg.bar2$any.promoter[df.gg.bar2$anchor1.promoter==T | df.gg.bar2$anchor2.promoter==T]="promoter"
df.promoter.perc=data.frame(group=rep(c("promoter","no-promoter"),3),frac=NA,group2=c(rep("common",2),rep(C1.lab,2),rep(C2.lab,2)),stringsAsFactors = F)
for (i in 1:nrow(df.promoter.perc)){
  df.promoter.perc$frac[i]=round(sum(df.gg.bar2$any.promoter==df.promoter.perc$group[i] & df.gg.bar2$group2==df.promoter.perc$group2[i])/sum(df.gg.bar2$group2==df.promoter.perc$group2[i]),digits = 2)
}
#write.table(df.promoter.perc,paste0(outdir,"/metrics_promoter_annot.tsv"),row.names = F,col.names = T,quote = F,sep="\t")

df.common$any.promoter="no-promoter"
df.common$any.promoter[df.common$anchor1.promoter==T | df.common$anchor2.promoter==T]="promoter"
df.promoter.perc.common=data.frame(group=rep(c("promoter","no-promoter"),3),frac=NA,group2=c(rep("stable",2),rep("increased",2),rep("decreased",2)),stringsAsFactors = F)
for (i in 1:nrow(df.promoter.perc.common)){
  df.promoter.perc.common$frac[i]=round(sum(df.common$any.promoter==df.promoter.perc.common$group[i] & df.common$loop.strength==df.promoter.perc.common$group2[i])/sum(df.common$loop.strength==df.promoter.perc.common$group2[i]),digits = 2)
}

# customize colors
df.gg2=df.gg
df.gg2$group2=as.character(df.gg2$group2)
df.gg2$group2[grep("common",df.gg2$group2)]="common"
df.gg2$group2=as.factor(df.gg2$group2)
df.gg2$group2=factor(df.gg2$group2,levels(df.gg2$group2))
ins.stat2=unique(as.character(df.gg2$group2))
df.gg2$group2=factor(df.gg2$group2,c("common",ins.stat2[ins.stat2 != "common"]))

keys.df3=data.frame(keys=levels(df.gg2$group2),color=NA)
keys.df3$color[keys.df3$keys =="common"]="darkgrey"
keys.df3$color[keys.df3$keys == C1.lab]="blue"
keys.df3$color[keys.df3$keys == C2.lab]="red"

df.gg$group2=as.factor(df.gg$group2)
df.gg$group2=factor(df.gg$group2,levels(df.gg$group2))
ins.stat2=unique(as.character(df.gg$group2))
df.gg$group2=factor(df.gg$group2,c(paste0("common-",C1.lab),paste0("common-",C2.lab),ins.stat2[grep("common",ins.stat2,invert = T)]))

keys.df2=data.frame(keys=levels(df.gg$group2),color=NA)
keys.df2$color[keys.df2$keys ==paste0("common-",C1.lab)]="darkgreen"
keys.df2$color[keys.df2$keys ==paste0("common-",C2.lab)]="black"
keys.df2$color[keys.df2$keys == C1.lab]="blue"
keys.df2$color[keys.df2$keys == C2.lab]="red"

df.gg$group=as.factor(df.gg$group)
df.gg$group=factor(df.gg$group,levels(df.gg$group))
ins.stat=unique(as.character(df.gg$group))
df.gg$group=factor(df.gg$group,c(ins.stat))

keys.df=data.frame(keys=levels(df.gg$group),color=NA)
keys.df$color[keys.df$keys == C1.lab]="blue"
keys.df$color[keys.df$keys == C2.lab]="red"

df.gg.common$group=as.factor(df.gg.common$group)
df.gg.common$group=factor(df.gg.common$group,levels(df.gg.common$group))
ins.stat=unique(as.character(df.gg.common$group))
df.gg.common$group=factor(df.gg.common$group,c(ins.stat))

keys.df.common=data.frame(keys=levels(df.gg.common$group),color=NA)
keys.df.common$color[keys.df.common$keys == C1.lab]="blue"
keys.df.common$color[keys.df.common$keys == C2.lab]="red"

df.common$loop.strength=as.factor(df.common$loop.strength)
df.common$loop.strength=factor(df.common$loop.strength,levels(df.common$loop.strength))
ins.stat=unique(as.character(df.common$loop.strength))
df.common$loop.strength=factor(df.common$loop.strength,c("stable",ins.stat[ins.stat != "stable"]))
df.common=df.common[order(df.common$loop.strength),]

keys.df.common2=data.frame(keys=levels(df.common$loop.strength),color=NA)
keys.df.common2$color[keys.df.common2$keys == "stable"]="darkgrey"
keys.df.common2$color[keys.df.common2$keys == "decreased"]="blue"
keys.df.common2$color[keys.df.common2$keys == "increased"]="red"

### Save plots in objects ###
p1=contactsByDist(df.gg,"group",is.common = F)
p2=contactsByDist(df.gg,"group2",is.common = F)
p3=contactsByDistByChr(df.gg,"group",is.common = F)
p4=contactsByDistByChr(df.gg,"group2",is.common = F)
p5=contactsByDist(df.gg.common,"group",is.common = T)
p6=contactsByDistByChr(df.gg.common,"group",is.common = T)
p7=lfcRankedContacts(df.common)
p8=hexbinContacts(df.common)
p9=scatterContactCount(df.common)
p10=distanceBox(df.common,max.dist=10000000)
p11=loopCountByDist(df.distance.num)
p12=countBarLoopClass(df.gg.bar2)
p13=countBarStacked(df.gg2)
p14=fractionBarStacked(df.perc)
p15=fractionBarStackedLoopClass(df.count)
p16=promoterCount(df.gg.bar2)
p17=promoterFraction(df.promoter.perc)
p18=promoterCountCommon(df.common)
p19=promoterFractionCommon(df.promoter.perc.common)

## ContactCounts ##
# loop contacts by distance
pdf(paste0(outdir,"/plots/loops_contactCounts_by_distance.pdf"), width=8 ,height=4) 
grid.arrange(arrangeGrob(p1,p2,nrow=1))
dev.off()

# loop contacts by distance by chromosome
pdf(paste0(outdir,"/plots/loops_contactCounts_by_distance_chrs.pdf"), width=10 ,height=7) 
p3
dev.off()

pdf(paste0(outdir,"/plots/loops_contactCounts_by_distance_loopClass_chrs.pdf"), width=10 ,height=7) 
p4
dev.off()

# common loops by distance / ranked (C1) vs lfc contactCounts
pdf(paste0(outdir,"/plots/common.loops_distance_contactCounts_lfc_ranked.pdf"), width=9 ,height=7) 
grid.arrange(arrangeGrob(p5,p7,nrow=1))
dev.off()

# common loop contacts by distance by chromosome
pdf(paste0(outdir,"/plots/common.loops_contactCounts_by_distance_chrs.pdf"), width=9 ,height=7) 
p6
dev.off()

# common loops scatter / hexbin contactCounts
pdf(paste0(outdir,"/plots/common.loops_contactCounts_scatter_hexbin.pdf"), width=8 ,height=4) 
grid.arrange(arrangeGrob(p9,p8,nrow=1))
dev.off()

## Distance ##
# distance boxplot / loopCount
pdf(paste0(outdir,"/plots/loop_distance_boxplot_loopCounts.pdf"), width=10 ,height=5) 
grid.arrange(arrangeGrob(p10,p11,nrow=1))
dev.off()

# barplots
pdf(paste0(outdir,"/plots/counts_perc_barplots.pdf"), width=10 ,height=7) 
grid.arrange(arrangeGrob(p12,p13,nrow=1),arrangeGrob(p15,p14,nrow=1),top="Common and group-specific loops")
dev.off()

# barplots promoter annotation
pdf(paste0(outdir,"/plots/counts_perc_promoterAnnot_barplots.pdf"), width=10 ,height=6) 
grid.arrange(arrangeGrob(p16,p17,nrow=1),arrangeGrob(p18,p19,nrow=1),top="Promoter-annotation: Common and group-specific loops",bottom="Criteria: one or more TSS in at least 1 of the loop-anchors")
dev.off()

#### Report ###
pdf(paste0(outdir,"/loops_report.pdf"), width=12 ,height=6) 
grid.arrange(arrangeGrob(p12,p13,nrow=1),arrangeGrob(p15,p14,nrow=1),top="Common and group-specific loops")
grid.arrange(arrangeGrob(p1,p2,nrow=1),bottom="Note: CPM normalized contactCounts were used.")
grid.arrange(arrangeGrob(p9,p8,nrow=1))
grid.arrange(arrangeGrob(p5,p7,nrow=1))
grid.arrange(arrangeGrob(p10,p11,nrow=1))
p3
#p4
p6
grid.arrange(arrangeGrob(p16,p17,nrow=1),arrangeGrob(p18,p19,nrow=1),top="Promoter-annotation: Common and group-specific loops",bottom="Criteria: one or more TSS in at least 1 of the loop-anchors")
dev.off()
