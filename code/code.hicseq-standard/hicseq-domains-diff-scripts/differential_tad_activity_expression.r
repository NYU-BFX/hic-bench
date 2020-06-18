#R
library(dplyr)
library(ggplot2)
library(GenomicRanges)

args <- commandArgs()
tad_activity_file <- args[3]
gene_tss_file <- args[4]
exprs_file_case_1 <- args[5]
exprs_file_case_2 <- args[6]
case_right <- args[7] #object1 -> S2
case_left <- args[8] #object2 -> S1 -> fold-change is S2/S1
bin.size <-  as.numeric(args[9])
min.TAD.size <- as.numeric(args[10])
logFC_threshold <- as.numeric(args[11])
fdr_threshold <- as.numeric(args[12])
meanDiff_threshold <- as.numeric(args[13])
out_prefix <- args[14]

# In communicating the tad_activity data.frame to text output, the scientific notation does not help.
options(scipen=999)

# cutoffs
tad_activity=read.csv(file=tad_activity_file,header=TRUE,sep="\t")
names(tad_activity)[6:7]=c("sample_2_mean","sample_1_mean") 
write.table(tad_activity,paste0(out_prefix,".tsv"),col.names = T,row.names = F,quote = F,sep="\t")

tad_activity <- tad_activity[(tad_activity$j - tad_activity$i) >= (min.TAD.size/bin.size),] # min-TAD size requirement
gene_tss <- read.csv(file=gene_tss_file, sep="\t", header=F)
gene_tss=gene_tss[,1:4]
names(gene_tss)=c("chr","start","end","geneID")

tad_activity$mean_overall <- apply(tad_activity[,c("sample_1_mean","sample_2_mean")],1,mean)
tad_activity_up <- (tad_activity[tad_activity$FDR <= fdr_threshold & tad_activity$logFC >= logFC_threshold & abs(tad_activity$mean_diff) >= meanDiff_threshold,])
tad_activity_down <- (tad_activity[tad_activity$FDR <= fdr_threshold & tad_activity$logFC <= -logFC_threshold & abs(tad_activity$mean_diff) >= meanDiff_threshold,])
tad_activity_unchanged <- (tad_activity[abs(tad_activity$mean_diff) < meanDiff_threshold | abs(tad_activity$logFC) < logFC_threshold |  tad_activity$FDR > fdr_threshold,])
tad_activity_unchanged_50 <- quantile(tad_activity$mean_overall, prob=0.5,na.rm=TRUE)
tad_activity_unchanged_low <- (tad_activity[abs(tad_activity$mean_diff) < meanDiff_threshold | abs(tad_activity$logFC) < logFC_threshold |  tad_activity$FDR > fdr_threshold &
                                              tad_activity$mean_overall < tad_activity_unchanged_50,])
tad_activity_unchanged_high <- tad_activity[abs(tad_activity$mean_diff) < meanDiff_threshold | abs(tad_activity$logFC) < logFC_threshold  |  tad_activity$FDR > fdr_threshold &
                                              tad_activity$mean_overall >= tad_activity_unchanged_50,]

write.table(file=paste(out_prefix,"_active-TADs.bed",sep=""),tad_activity_up[,c(1,3,5)],sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(file=paste(out_prefix,"_inactive-TADs.bed",sep=""),tad_activity_down[,c(1,3,5)],sep="\t",quote=FALSE,row.names=FALSE, col.names=FALSE)
write.table(file=paste(out_prefix,"_unchanged-TADs.bed",sep=""),tad_activity_unchanged[,c(1,3,5)],sep="\t",quote=FALSE,row.names=FALSE, col.names=FALSE)
write.table(file=paste(out_prefix,"_active-TADs.tsv",sep=""),tad_activity_up,sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(file=paste(out_prefix,"_inactive-TADs.tsv",sep=""),tad_activity_down,sep="\t",quote=FALSE,row.names=FALSE, col.names=FALSE)
write.table(file=paste(out_prefix,"_unchanged-TADs.tsv",sep=""),tad_activity_unchanged,sep="\t",quote=FALSE,row.names=FALSE, col.names=FALSE)

# classify tads activity-wise
tad_activity$activity.group[tad_activity$FDR <= fdr_threshold & tad_activity$logFC >= logFC_threshold & abs(tad_activity$mean_diff) >= 0.1]="increased activity"
tad_activity$activity.group[tad_activity$FDR <= fdr_threshold & tad_activity$logFC <= -logFC_threshold & abs(tad_activity$mean_diff) >= 0.1]="decreased activity"
tad_activity$activity.group[abs(tad_activity$logFC) < logFC_threshold | tad_activity$FDR > fdr_threshold | abs(tad_activity$mean_diff) < 0.1]="stable TADs"
tad_activity$activity.group=as.factor(tad_activity$activity.group)
tad_activity$activity.group=factor(tad_activity$activity.group,levels(tad_activity$activity.group))
ins.stat=unique(as.character(tad_activity$activity.group))
tad_activity$activity.group=factor(tad_activity$activity.group,c("stable TADs",ins.stat[ins.stat != "stable TADs"]))
tad_activity=tad_activity[order(tad_activity$activity.group),]

# annotate gene-TADs
tad_activity$geneIDs=NA
tad_activity$start=tad_activity$TAD_start
tad_activity$end=tad_activity$TAD_end

for (tad.num in 1:nrow(tad_activity)){ # extend TAD coordinates to cover the boundaries
  if(tad.num > 1 & tad.num < nrow(tad_activity)){
    if(tad_activity$chr[tad.num-1] == tad_activity$chr[tad.num]) { tad_activity$start[tad.num]=tad_activity$TAD_end[tad.num-1] }
    if(tad_activity$chr[tad.num+1] == tad_activity$chr[tad.num]) { tad_activity$end[tad.num]=tad_activity$TAD_start[tad.num+1] }
  } else {
    if (tad.num == 1){ 
      if(tad_activity$chr[tad.num+1] == tad_activity$chr[tad.num]) { tad_activity$end[tad.num]=tad_activity$TAD_start[tad.num+1] }
    }
    if (tad.num == nrow(tad_activity)){ 
      if(tad_activity$chr[tad.num-1] == tad_activity$chr[tad.num]) { tad_activity$start[tad.num]=tad_activity$TAD_end[tad.num-1] }
    }
  }
}

gene_tss.gr=makeGRangesFromDataFrame(gene_tss,keep.extra.columns = T)
tad_activity.gr=makeGRangesFromDataFrame(tad_activity,keep.extra.columns = T)

for (tad.num in 1:nrow(tad_activity)) {
  gene_hits=subsetByOverlaps(gene_tss.gr,tad_activity.gr[tad.num])
  tad_activity$geneIDs[tad.num]=paste(gene_hits$geneID,collapse = ",")
}
write.table(file=paste(out_prefix,"_annotated.tsv",sep=""),tad_activity,sep="\t",quote=FALSE,row.names=FALSE)

# Same output in .bed format, with the fold change as 5th column
# Notice the length shortening by 500 bp both sides. This is only to ensure that pygenometracks draws tracks in a single line.
# Do not directly use this final_results_forpygenometracks.bed for downstream analysis. For visualization on pygenometracks only.
write.table(data.frame(tad_activity$chr, tad_activity$TAD_start + 500, tad_activity$TAD_end - 500, "Domains_diff_FC", ifelse(!is.na(tad_activity$logFC), tad_activity$logFC, 0.0), "*"), file=paste0(out_prefix, "_forpygenometracks.bed"), sep="\t", row.names=F, col.names=F, quote=F)

# get max absolute lfc & fdr values (helps to set axes limits)
tad_activity.clean=tad_activity[abs(tad_activity$mean_diff) >= meanDiff_threshold &
                                  is.finite(tad_activity$sample_1_mean) &
                                  is.finite(tad_activity$sample_2_mean) &
                                  is.finite(tad_activity$logFC) &
                                  is.finite(-log(tad_activity$FDR)),]

abs.max.x=max(abs(tad_activity.clean$logFC),na.rm = T)
abs.max.y=max(abs(-log(tad_activity.clean$FDR)),na.rm = T)
box.max.y=max(abs(c(tad_activity.clean$sample_1_mean,tad_activity.clean$sample_2_mean)),na.rm = T)

# get label character count to adjust margins
bottom.mar=5 #default bottom margin 
right.mar=15 #default right margin
max.char=max(c(nchar(case_right),nchar(case_left)),na.rm = T) #max label character value
char.diff.r=max.char-9 #9 is the max num of label characters that fit using the default right margin value
char.diff.b=max.char-7 #7 is the max num of label characters that fit using the default bottom margin value

mcf=0.55 # margin correction factor: total extra characters (char.diff) * mcf
if (char.diff.b > 0){b=bottom.mar+(char.diff.b*mcf)} else {b=bottom.mar} 
if (char.diff.r > 0){r=right.mar+(char.diff.r*mcf)} else {r=right.mar}

### Mean TAD activity BOXPLOTS ###
png(file=paste(out_prefix, "_mean-TAD_activity.png", sep=""),width=3072, height=2048, pointsize=60)
par(mar=c(b,4,1.1,r),xpd=T)
boxplot(tad_activity_up$sample_2_mean,
        tad_activity_up$sample_1_mean,
        tad_activity_unchanged_low$sample_2_mean,
        tad_activity_unchanged_low$sample_1_mean,
        tad_activity_unchanged_high$sample_2_mean,
        tad_activity_unchanged_high$sample_1_mean,
        tad_activity_down$sample_1_mean,
        tad_activity_down$sample_2_mean,
        beside=TRUE, col=c("darkred","darkred","grey","grey","grey","grey","blue","blue"),bty="n",axes=FALSE,xlab="",ylab="Mean TAD activity",
        lwd=7,pch=20)
box(lwd=7)
axis(side=1,lwd=7,labels=c(case_right,case_left,case_right,case_left,case_right,case_left,case_left,case_right),at=seq(1,8,1),las=2)
axis(side=2,lwd=7,las=2)
legend(x=9,y=box.max.y*0.7,legend = c(paste0("increased activity in ",case_right),paste0("decreased activity in ",case_right),"stable TADs (Low)","stable TADs (High)"),
       pch=22,col=c("black"),pt.bg=c("darkred","blue","grey","grey"),cex=1,box.col = NA,pt.cex = 2,pt.lwd = 4)
dev.off()

### VOLCANO PLOT ###
#adjust margins
right.mar=10 #default right margin
mcf=0.4 
char.diff.r=max.char-8 #8 is the max num of label characters that fit using the default right margin value
if (char.diff.r > 0){r=right.mar+(char.diff.r*mcf)} else {r=right.mar}

png(file=paste(out_prefix, "_volcano_TADs.png", sep=""),width=3755, height=2048, pointsize=110)
par(mar=c(4.2,4,2,r),xpd=T)
plot(tad_activity$logFC, -log(tad_activity$FDR), xlab="Mean TAD activity logFC",ylab="-log(FDR)",pch=20,col="light grey",
     cex=1,bty="n",main = paste0(case_right, " vs ",case_left),cex.main=0.6,xlim = c(-abs.max.x,abs.max.x),ylim=c(0,abs.max.y*1.1))
box(lwd=8)
axis(side = 1, lwd = 8)
axis(side = 2, lwd = 8)
points(tad_activity_unchanged$logFC, -log(tad_activity_unchanged$FDR),col="grey",pch=20,cex=1)
points(tad_activity_up$logFC, -log(tad_activity_up$FDR),col="darkred",pch=20,cex=1)
points(tad_activity_down$logFC, -log(tad_activity_down$FDR),col="blue",pch=20,cex=1)
abline(h=-log(fdr_threshold), col="red", lwd=8, lty=2,xpd=F)
abline(v=-logFC_threshold,col="red",lwd=8,lty=2,xpd=F)
abline(v=logFC_threshold,col="red",lwd=8,lty=2,xpd=F)
legend(x=abs.max.x*1.15,y=abs.max.y*0.7,legend = c(paste0("increased activity in ",case_right),paste0("decreased activity in ",case_right),("stable TADs")),
       pch=21,col=c("black"),pt.bg=c("darkred","blue","grey"),cex=0.7,box.col = NA,pt.cex = 1,pt.lwd = 4)
dev.off()

### SCATTER PLOT ###
keys.df=data.frame(keys=levels(tad_activity$activity.group),color=NA)
keys.df$color[keys.df$keys =="stable TADs"]="grey"
keys.df$color[keys.df$keys =="increased activity"]="darkred"
keys.df$color[keys.df$keys =="decreased activity"]="blue"

pdf(file=paste0(out_prefix,"_scatter_meanTadActivity.pdf"),width = 10,height = 6)
ggplot(tad_activity, aes(sample_1_mean, sample_2_mean,color=activity.group))+ 
  geom_point(size=3)+
  scale_color_manual(values=keys.df$color)+
  geom_abline(slope = 1,intercept = 0,linetype=3,color="black",size=1)+
  xlab(paste0("Mean TAD activity [",case_left,"]"))+
  ylab(paste0("Mean TAD activity [",case_right,"]"))+
  theme_classic(base_size=14)+
  theme(legend.title=element_blank())
dev.off()
