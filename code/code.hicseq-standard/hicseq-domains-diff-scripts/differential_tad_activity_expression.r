#R
library(dplyr)

args <- commandArgs()

tad_activity_file <- args[3]
gene_bins_file <- args[4]
exprs_file_case_1 <- args[5]
exprs_file_case_2 <- args[6]
case_left <- args[7]
case_right <- args[8]
bin.size <-  as.numeric(args[9])
min.TAD.size <- as.numeric(args[10])
out_prefix <- args[11]

tad_activity <- read.csv(file=tad_activity_file,header=TRUE,sep="\t")
# min-TAD size requirement
tad_activity <- tad_activity[(tad_activity$j - tad_activity$i) >= (min.TAD.size/bin.size),]
#tad_activity <- tad_activity[tad_activity$logFC != Inf & tad_activity$logFC != -Inf ,]
gene_bins <- read.csv(file=gene_bins_file, sep="\t", header=TRUE)

# necessary??
gene_bins <- gene_bins[gene_bins$type=="protein_coding" | gene_bins$type=="processed_transcript",]

logFC_threshold <- 0.3
fdr_threshold <- 0.1

tad_activity$mean_overall <- apply(tad_activity[,c("sample_1_mean","sample_2_mean")],1,mean)
tad_activity_up <- (tad_activity[tad_activity$FDR <= fdr_threshold & tad_activity$mean_diff >= logFC_threshold,])
tad_activity_down <- (tad_activity[tad_activity$FDR <= fdr_threshold & tad_activity$mean_diff <= -logFC_threshold,])
#tad_activity_unchanged <- (tad_activity[tad_activity$FDR >= fdr_threshold | abs(tad_activity$mean_diff) < logFC_threshold ,])
tad_activity_unchanged <- (tad_activity[tad_activity$mean_diff >= -0.1 & (tad_activity$mean_diff) <= 0.1,])
tad_activity_unchanged_50 <- quantile(tad_activity$mean_overall, prob=0.5,na.rm=TRUE)
tad_activity_unchanged_low <- tad_activity[tad_activity$mean_diff >= -0.1 & (tad_activity$mean_diff) <= 0.1 &
        tad_activity$mean_overall < tad_activity_unchanged_50,]
tad_activity_unchanged_high <- tad_activity[tad_activity$mean_diff >= -0.1 & (tad_activity$mean_diff) <= 0.1 &
        tad_activity$mean_overall >= tad_activity_unchanged_50,]

write.table(file=paste(out_prefix,"_active-TADs.bed",sep=""),tad_activity_up[,c(1,3,5)],sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(file=paste(out_prefix,"_inactive-TADs.bed",sep=""),tad_activity_down[,c(1,3,5)],sep="\t",quote=FALSE,row.names=FALSE, col.names=FALSE)
write.table(file=paste(out_prefix,"_unchanged-TADs.bed",sep=""),tad_activity_unchanged[,c(1,3,5)],sep="\t",quote=FALSE,row.names=FALSE, col.names=FALSE)

write.table(file=paste(out_prefix,"_active-TADs.tsv",sep=""),tad_activity_up,sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(file=paste(out_prefix,"_inactive-TADs.tsv",sep=""),tad_activity_down,sep="\t",quote=FALSE,row.names=FALSE, col.names=FALSE)
write.table(file=paste(out_prefix,"_unchanged-TADs.tsv",sep=""),tad_activity_unchanged,sep="\t",quote=FALSE,row.names=FALSE, col.names=FALSE)

png(file=paste(out_prefix, "_mean-TAD_activity.png", sep=""),width=2048, height=2048, pointsize=60)
par(mar=c(4,4,0.1,0.1))
boxplot(tad_activity_up$sample_1_mean,
        tad_activity_up$sample_2_mean,
        tad_activity_down$sample_1_mean,
        tad_activity_down$sample_2_mean,
        tad_activity_unchanged_low$sample_1_mean,
        tad_activity_unchanged_low$sample_2_mean,
        tad_activity_unchanged_high$sample_1_mean,
        tad_activity_unchanged_high$sample_2_mean,
        beside=TRUE, col=c("darkred","darkred","blue","blue","grey","grey","grey","grey"),bty="n",axes=FALSE,xlab="",ylab="Mean TAD activity",
        lwd=7,pch=20)
box(lwd=7)
axis(side=1,lwd=7,labels=c("Sample 1","Sample 2","Sample 1","Sample 2","Sample 1","Sample 2","Sample 1","Sample 2"),at=seq(1,8,1),las=2)
axis(side=2,lwd=7,las=2)
dev.off()

png(file=paste(out_prefix, "_volcano_TADs.png", sep=""),width=2048, height=2048, pointsize=110)
#max_x <- max(abs(tad_activity$logFC),na.rm=TRUE)
par(mar=c(4.2,4,0.4,0.4))
plot(tad_activity$mean_diff, -log(tad_activity$FDR), xlab="log(fold-change)",ylab="-log(FDR)",pch=20,col="light grey",
	cex=1.5,bty="n")
#XX
#max_x <- 3.5
#max_y <- 175
#plot(tad_activity$logFC, -log(tad_activity$FDR), xlab="log(fold-change)",ylab="-log(FDR)",pch=20,col="grey",
#	cex=1.5,xlim=c(((-1)*max_x),max_x), ylim=c(0,max_y),bty="n")

#font=2,font.lab=2,
box(lwd=8)
axis(side = 1, lwd = 8)
axis(side = 2, lwd = 8)
points(tad_activity_up$mean_diff, -log(tad_activity_up$FDR),col="darkred",pch=20,cex=1.5)
points(tad_activity_down$mean_diff, -log(tad_activity_down$FDR),col="blue",pch=20,cex=1.5)
# Added unchanging TADs too
points(tad_activity_unchanged$mean_diff, -log(tad_activity_unchanged$FDR),col="grey",pch=20,cex=1.5)
abline(h=-log2(fdr_threshold), col="red", lwd=8, lty=2)
abline(v=-logFC_threshold,col="red",lwd=8,lty=2)
abline(v=logFC_threshold,col="red",lwd=8,lty=2)
dev.off()

png(file=paste(out_prefix, "_volcano_TADs_legend.png", sep=""),width=2048, height=2048, pointsize=110)
plot.new()
legend("center",legend=c("Higher activity in T-ALL", "Less activity in T-ALL"),fill=c("darkred","blue"),lwd=8)
dev.off()

#tad_extension <- 0
#tad_extension <- 120000
tad_extension <- 40000

nrow(tad_activity)

#gene_bin_hits_exprs_up <- data.frame()
gene_hits <- vector()
for (i in 1:nrow(tad_activity)) {
	gene_bin_hits <- gene_bins[as.character(gene_bins$chr) == as.character(tad_activity[i,]$chr) & gene_bins$start >= (tad_activity[i,]$TAD_start - tad_extension) & 
		gene_bins$end <= (tad_activity[i,]$TAD_end + tad_extension),]
	if (nrow(gene_bin_hits)==0) {
		gene_hits <- rbind(gene_hits,"NA")
	} else {
		gene_hits <- rbind(gene_hits,paste(as.character(gene_bin_hits$symbol),collapse=","))
	}
#		gene_hits_up <- c(gene_hits_up, as.character(gene_bin_hits$ID))
#		gene_bin_hits_exprs <- rbind(gene_bin_hits_exprs, exprs[rownames(exprs) %in% gene_bin_hits$ID,])
}
tad_activity <- cbind(tad_activity,gene_hits)
write.table(file=paste(out_prefix,"_annotated.tsv",sep=""),tad_activity,sep="\t",quote=FALSE,row.names=FALSE)

if (exprs_file_case_1 == "FALSE") {
	stop(0)
}

# Read gene expression data and filter for expressed genes
exprs_case_1 <- as.data.frame(readRDS(exprs_file_case_1))
rownames(exprs_case_1) <- as.character(exprs_case_1$name)
exprs_case_2 <- as.data.frame(readRDS(exprs_file_case_2))
rownames(exprs_case_2) <- as.character(exprs_case_2$name)
filtered_names <- exprs_case_1$cpm > 1 | exprs_case_2$cpm > 1

exprs_case_1 <- exprs_case_1[filtered_names,]
exprs_case_2 <- exprs_case_2[filtered_names,]

gene_bin_hits_exprs_up <- data.frame()
gene_hits_up <- vector()
for (i in 1:nrow(tad_activity_up)) {
		gene_bin_hits <- gene_bins[as.character(gene_bins$chr) == as.character(tad_activity_up[i,]$chr) & gene_bins$start >= (tad_activity_up[i,]$TAD_start - tad_extension) & 
			gene_bins$end <= (tad_activity_up[i,]$TAD_end + tad_extension),]
		gene_hits_up <- c(gene_hits_up, as.character(gene_bin_hits$ID))
#		gene_bin_hits_exprs_up <- rbind(gene_bin_hits_exprs_up, exprs[rownames(exprs) %in% gene_bin_hits$ID,])
		
}

gene_bin_hits_exprs_down <- data.frame()
gene_hits_down <- vector()
for (i in 1:nrow(tad_activity_down)) {
		gene_bin_hits <- gene_bins[as.character(gene_bins$chr) == as.character(tad_activity_down[i,]$chr) & gene_bins$start >= (tad_activity_down[i,]$TAD_start - tad_extension) & 
			gene_bins$end <= (tad_activity_down[i,]$TAD_end + tad_extension),]
		gene_hits_down <- c(gene_hits_down, as.character(gene_bin_hits$ID))
#		gene_bin_hits_exprs_down <- rbind(gene_bin_hits_exprs_down, as.data.frame(exprs[rownames(exprs) %in% gene_bin_hits$ID,]))
}

# Added unchanging TADs here too
gene_bin_hits_exprs_unchanged <- data.frame()
gene_hits_unchanged <- vector()
for (i in 1:nrow(tad_activity_unchanged)) {
		gene_bin_hits <- gene_bins[as.character(gene_bins$chr) == as.character(tad_activity_unchanged[i,]$chr) & gene_bins$start >= (tad_activity_unchanged[i,]$TAD_start - tad_extension) & 
			gene_bins$end <= (tad_activity_unchanged[i,]$TAD_end + tad_extension),]
		gene_hits_unchanged <- c(gene_hits_unchanged, as.character(gene_bin_hits$ID))
#		gene_bin_hits_exprs_up <- rbind(gene_bin_hits_exprs_up, exprs[rownames(exprs) %in% gene_bin_hits$ID,])
		
}
gene_hits_up <- na.omit(gene_hits_up)
gene_hits_down <- na.omit(gene_hits_down)
gene_hits_unchanged <- na.omit(gene_hits_unchanged)
#nrow(gene_hit_down)
#ERG
#tad_activity[tad_activity$chr=="chr21" & tad_activity$TAD_start < 40000001 & tad_activity$TAD_end > 40040000,]

write.table(file=paste(out_prefix,"_genes-in-active-TADs.tsv",sep=""),(gene_hits_up),sep="\t",quote=FALSE,row.names=FALSE)
write.table(file=paste(out_prefix,"_genes-in-inactive-TADs.tsv",sep=""),(gene_hits_down),sep="\t",quote=FALSE,row.names=FALSE)
write.table(file=paste(out_prefix,"_genes-in-unchanged-TADs.tsv",sep=""),(gene_hits_unchanged),sep="\t",quote=FALSE,row.names=FALSE)

#rpkm_cd4_up <- data.frame(rpkm_cd4_up=apply(gene_bin_hits_exprs_up[,4:16],1,mean), name=rownames(gene_bin_hits_exprs_up))
#rpkm_cd4_down <- data.frame(rpkm_cd4_down=apply(gene_bin_hits_exprs_down[,4:16],1,mean), name=rownames(gene_bin_hits_exprs_down))

rpkm_left_up <- exprs_case_1[rownames(exprs_case_1) %in% gene_hits_up,]
rpkm_left_down <- exprs_case_1[exprs_case_1$name %in% gene_hits_down,]
rpkm_right_up <- exprs_case_2[exprs_case_2$name %in% gene_hits_up,]
rpkm_right_down <- exprs_case_2[exprs_case_2$name %in% gene_hits_down,]
rpkm_left_unchanged <- exprs_case_1[exprs_case_1$name %in% gene_hits_unchanged,]
rpkm_right_unchanged <- exprs_case_2[exprs_case_2$name %in% gene_hits_unchanged,]

# left = sample1
# right = sample2

str(rpkm_left_up)
#rpkm_cd4_up <- data.frame(rpkm_cd4_up=apply(gene_bin_hits_exprs_up$rpkm_mean_cd4,1,mean), name=rownames(gene_bin_hits_exprs_up))
#rpkm_cd4_down <- data.frame(rpkm_cd4_down=apply(gene_bin_hits_exprs_down$rpkm_mean_cd4,1,mean), name=rownames(gene_bin_hits_exprs_down))
#rpkm_tall_up <- data.frame(rpkm_tall_up=gene_bin_hits_exprs_up[,case_left], name=rownames(gene_bin_hits_exprs_up))
#rpkm_tall_down <- data.frame(rpkm_tall_down=gene_bin_hits_exprs_down[,case_left], name=rownames(gene_bin_hits_exprs_down))

# I changed rpkm with cpm

# Add pseudocount at logFC calculations
psd <- 0.1
rpkm_merged_up <- merge(rpkm_left_up, rpkm_right_up, by="name",all=FALSE)
#logFC_up <- log2(rpkm_merged_up$cpm.x / rpkm_merged_up$cpm.y)
logFC_up <- log2((rpkm_merged_up$cpm.x +psd)/ (rpkm_merged_up$cpm.y + psd))
rpkm_merged_down <- merge(rpkm_left_down, rpkm_right_down, by="name",all=FALSE)
#logFC_down <- log2(rpkm_merged_down$cpm.x / rpkm_merged_down$cpm.y)
logFC_down <- log2((rpkm_merged_down$cpm.x +psd)/ (rpkm_merged_down$cpm.y+psd))
rpkm_merged_unchanged <- merge(rpkm_left_unchanged, rpkm_right_unchanged, by="name",all=FALSE)
#logFC_unchanged <- log2(rpkm_merged_unchanged$cpm.x / rpkm_merged_unchanged$cpm.y)
logFC_unchanged <- log2((rpkm_merged_unchanged$cpm.x +psd)/ (rpkm_merged_unchanged$cpm.y+psd))

#logFC_down

# Create a data frame with the logFC expression of all common genes between the samples
rpkm_merged_all_genes <- merge(exprs_case_1, exprs_case_2, by = "name", all=FALSE)
#logFC_all_genes <- log2(rpkm_merged_all_genes$cpm.x / rpkm_merged_all_genes$cpm.y)
logFC_all_genes <- log2((rpkm_merged_all_genes$cpm.x +psd)/ (rpkm_merged_all_genes$cpm.y+psd))

# Create data frames with the genes that have logFC > 0.5 for the active TADs
#			and print them

rpkm_names_merged_up <- as.data.frame(rpkm_merged_up[,c("name")])
colnames(rpkm_names_merged_up)[1] <- "names"
rpkm_names_merged_up$logFC <- logFC_up
rpkm_names_merged_up <- rpkm_names_merged_up[rpkm_names_merged_up$logFC > 0.5,]

rpkm_names_merged_down <- as.data.frame(rpkm_merged_down[,c("name")])
colnames(rpkm_names_merged_down)[1] <- "names"
rpkm_names_merged_down$logFC <- logFC_down
rpkm_names_merged_down <- rpkm_names_merged_down[rpkm_names_merged_down$logFC < -0.5,]

rpkm_names_merged_unchanged <- as.data.frame(rpkm_merged_unchanged[,c("name")])
colnames(rpkm_names_merged_unchanged)[1] <- "names"
rpkm_names_merged_unchanged$logFC <- logFC_unchanged
rpkm_names_merged_unchanged <- rpkm_names_merged_unchanged[abs(rpkm_names_merged_unchanged$logFC)<0.5,]

write.table(file=paste(out_prefix,"_genes-in-active-TADs_logFC.tsv",sep=""),rpkm_names_merged_up,sep="\t",quote=FALSE,row.names=FALSE)
write.table(file=paste(out_prefix,"_genes-in-inactive-TADs_logFC.tsv",sep=""),rpkm_names_merged_down,sep="\t",quote=FALSE,row.names=FALSE)
write.table(file=paste(out_prefix,"_genes-in-unchanged-TADs_logFC.tsv",sep=""),rpkm_names_merged_unchanged,sep="\t",quote=FALSE,row.names=FALSE)

#logFC_up <- log2(rpkm_right_up$cpm / rpkm_left_up$cpm)
#logFC_down <- log2(rpkm_right_down$cpm / rpkm_left_down$cpm)
#logFC_unchanged <- log2(rpkm_right_unchanged$cpm / rpkm_left_unchanged$cpm)


#str(logFC_up)

#rpkm_right_up$logFC  <- log2(rpkm_right_up$cpm / rpkm_left_up$cpm)
#rpkm_right_down$logFC <- log2(rpkm_right_down$cpm / rpkm_left_down$cpm)
#rpkm_right_unchanged$logFC <- log2(rpkm_right_unchanged$cpm / rpkm_left_unchanged$cpm)

#write.table(file=paste(out_prefix,"_genes-in-active-TADs-logFC.tsv",sep=""),cbind(as.character(rpkm_right_up$name),rpkm_right_up$logFC),sep="\t",quote=FALSE,row.names=FALSE)
#write.table(file=paste(out_prefix,"_genes-in-inactive-TADs-logFC.tsv",sep=""),cbind(as.character(rpkm_right_down$name),rpkm_right_down$logFC),sep="\t",quote=FALSE,row.names=FALSE)
#write.table(file=paste(out_prefix,"_genes-in-unchanged-TADs-logFC.tsv",sep=""),cbind(as.character(rpkm_right_unchanged$name),rpkm_right_unchanged$logFC),sep="\t",quote=FALSE,row.names=FALSE)

#logFC_up_cleaned <- logFC_up[logFC_up != Inf & logFC_up != -Inf]
logFC_up_cleaned <- logFC_up

if (length(logFC_up_cleaned) != 0){
logFC_up_cleaned_test <- t.test(logFC_up_cleaned, na.action="na.omit", alternative="greater")
#logFC_down_cleaned <- logFC_down[logFC_down != Inf & logFC_down != -Inf]
logFC_down_cleaned <- logFC_down
logFC_down_cleaned_test <- t.test(logFC_down_cleaned, na.action="na.omit", alternative="less")
#logFC_unchanged_cleaned <- logFC_unchanged[logFC_unchanged != Inf & logFC_unchanged != -Inf]
logFC_unchanged_cleaned <- logFC_unchanged
logFC_unchanged_cleaned_test <- t.test(logFC_unchanged_cleaned, na.action="na.omit", alternative="two.sided")

#logFC_all_genes_cleaned <- logFC_all_genes[logFC_all_genes != Inf & logFC_all_genes != -Inf]
logFC_all_genes_cleaned <- logFC_all_genes

logFC_up_unchanged_cleaned_test <- t.test(logFC_up_cleaned, logFC_unchanged_cleaned, na.action="na.omit", alternative="greater")
logFC_down_unchanged_cleaned_test <- t.test(logFC_down_cleaned, logFC_unchanged_cleaned, na.action="na.omit", alternative="less")
logFC_up_down_cleaned_test <- t.test(logFC_up_cleaned, logFC_down_cleaned, na.action="na.omit", alternative="greater")
logFC_all_genes_unchanged_cleaned_test <- t.test(logFC_all_genes_cleaned, logFC_unchanged_cleaned, na.action="na.omit", alternative = "two.sided")
 
# ADD OTHER p-values TOO and print them
test_results <- rbind(cbind(logFC_up_cleaned_test$estimate, logFC_up_cleaned_test$p.value),cbind(logFC_unchanged_cleaned_test$estimate, logFC_unchanged_cleaned_test$p.value), cbind(logFC_down_cleaned_test$estimate, logFC_down_cleaned_test$p.value),cbind(logFC_up_unchanged_cleaned_test$estimate, logFC_up_unchanged_cleaned_test$p.value), cbind(logFC_down_unchanged_cleaned_test$estimate, logFC_down_unchanged_cleaned_test$p.value), cbind(logFC_up_down_cleaned_test$estimate, logFC_up_down_cleaned_test$p.value), cbind(logFC_all_genes_unchanged_cleaned_test$estimate, logFC_all_genes_unchanged_cleaned_test$p.value))

rownames(test_results) <- c("Higher TAD activity", "Same TAD activity", "Less TAD activity", "Higher-Same TAD activity, High mean","Higher-Same TAD activity, Same mean", "Less-Same TAD activity, Less mean", "Less-Same TAD activity, Same mean","Higher-Less TAD activity, High mean", "Higher-Less TAD activity, Less mean", "All genes-unchanged, Mean All", "All genes-Unchanged, Mean unchanged")
colnames(test_results) <- c("Mean differences from 0", "p-value")
write.table(file=paste(out_prefix,"_logFC-differences_pvalues.tsv",sep=""),test_results,quote=FALSE, row.names=TRUE, sep="\t")

# compare up vs unchanged
# compare down vs unchanged 
# compare up vs down
# p values around 1 -> wrong -> change greater and less
# 2 sample tests

png(file=paste(out_prefix, "_boxplot_RNAexprs_logFC.png",sep=""),width=1200,height=2048,pointsize=110)
box_names <- c(case_right,case_left)
#par(mfrow=c(1,1),mar=c(8,4,4,4))
#y_max <- max(c(rpkm_left_up,gene_bin_hits_exprs_up[,tall_exprs_name],rpkm_left_down,gene_bin_hits_exprs_down[,tall_exprs_name]),
#	na.rm=TRUE)
#y_min <- min(c(rpkm_left_up,gene_bin_hits_exprs_up[,tall_exprs_name],rpkm_left_down,gene_bin_hits_exprs_down[,tall_exprs_name]),
#	na.rm=TRUE)
y_min <- -6
y_max <- 6
par(mar=c(0.1,4,0.1,0.1))

# Add the unchanged in he middle of down and up
# Also, add the logFC of all genes as an additional control
boxplot(logFC_down, logFC_unchanged, logFC_up, logFC_all_genes,beside=TRUE,  notch=TRUE,col=c("blue","grey","darkred","dimgrey"),las=2,
	ylab="RNA logFC",na.action="na.omit", ylim=c(y_min,y_max),byt="n",pch=20,font=2,lwd=8,outline=FALSE)
#boxplot(logFC_down, logFC_up, beside=TRUE, notch=TRUE,col=c("blue","darkred"),na.action="na.omit", add=TRUE,las=2, ylim=c(y_min,y_max),bty="n",pch=20,lwd=8)
abline(h=0, lwd=10,col="black")

stripchart(list(logFC_down=logFC_down,logFC_unchanged=logFC_unchanged, logFC_up=logFC_up, logFC_all_genes=logFC_all_genes),
        col=c("blue","grey","darkred","dimgrey"),add=TRUE,pch=20,cex=1.4,vertical=TRUE,method="jitter")

boxplot(logFC_down, logFC_unchanged, logFC_up, logFC_all_genes,beside=TRUE, notch=TRUE,col=c("blue","grey","darkred","dimgrey"),na.action="na.omit", add=TRUE,las=2, ylim=c(y_min,y_max),bty="n",pch=20,lwd=8,outline=FALSE)
box(lwd=8)
axis(side = 2, lwd = 8,labels=FALSE)
#legend("topright", legend=c("Less TAD activity","Higher TAD activity"), fill=c("blue","darkred"))
#boxplot(rpkm_left_down,gene_bin_hits_exprs_down[,tall_exprs_name],beside=TRUE, main="Less TAD activity in T-ALL",notch=TRUE,col="lightgray",las=2,
#	names=box_names,ylab="logFPKM",ylim=c(y_min,y_max))
dev.off()
} else {

logFC_down_cleaned <- logFC_down
logFC_down_cleaned_test <- t.test(logFC_down_cleaned, na.action="na.omit", alternative="less")
#logFC_unchanged_cleaned <- logFC_unchanged[logFC_unchanged != Inf & logFC_unchanged != -Inf]
logFC_unchanged_cleaned <- logFC_unchanged
logFC_unchanged_cleaned_test <- t.test(logFC_unchanged_cleaned, na.action="na.omit", alternative="two.sided")

#logFC_all_genes_cleaned <- logFC_all_genes[logFC_all_genes != Inf & logFC_all_genes != -Inf]
logFC_all_genes_cleaned <- logFC_all_genes

#logFC_up_unchanged_cleaned_test <- t.test(logFC_up_cleaned, logFC_unchanged_cleaned, na.action="na.omit", alternative="greater")
logFC_down_unchanged_cleaned_test <- t.test(logFC_down_cleaned, logFC_unchanged_cleaned, na.action="na.omit", alternative="less")
#logFC_up_down_cleaned_test <- t.test(logFC_up_cleaned, logFC_down_cleaned, na.action="na.omit", alternative="greater")
logFC_all_genes_unchanged_cleaned_test <- t.test(logFC_all_genes_cleaned, logFC_unchanged_cleaned, na.action="na.omit", alternative = "two.sided")

# ADD OTHER p-values TOO and print them
#test_results <- rbind(cbind(logFC_up_cleaned_test$estimate, logFC_up_cleaned_test$p.value),cbind(logFC_unchanged_cleaned_test$estimate, logFC_unchanged_cleaned_test$p.value), cbind(logFC_down_cleaned_test$estimate, logFC_down_cleaned_test$p.value),cbind(logFC_up_unchanged_cleaned_test$estimate, logFC_up_unchanged_cleaned_test$p.value), cbind(logFC_down_unchanged_cleaned_test$estimate, logFC_down_unchanged_cleaned_test$p.value), cbind(logFC_up_down_cleaned_test$estimate, logFC_up_down_cleaned_test$p.value), cbind(logFC_all_genes_unchanged_cleaned_test$estimate, logFC_all_genes_unchanged_cleaned_test$p.value))

#rownames(test_results) <- c("Higher TAD activity", "Same TAD activity", "Less TAD activity", "Higher-Same TAD activity, High mean","Higher-Same TAD activity, Same mean", "Less-Same TAD activity, Less mean", "Less-Same TAD activity, Same mean","Higher-Less TAD activity, High mean", "Higher-Less TAD activity, Less mean", "All genes-unchanged, Mean All", "All genes-Unchanged, Mean unchanged")
#colnames(test_results) <- c("Mean differences from 0", "p-value")
#write.table(file=paste(out_prefix,"_logFC-differences_pvalues.tsv",sep=""),test_results,quote=FALSE, row.names=TRUE, sep="\t")

png(file=paste(out_prefix, "_boxplot_RNAexprs_logFC.png",sep=""),width=1200,height=2048,pointsize=110)
box_names <- c(case_right,case_left)
#par(mfrow=c(1,1),mar=c(8,4,4,4))
#y_max <- max(c(rpkm_left_up,gene_bin_hits_exprs_up[,tall_exprs_name],rpkm_left_down,gene_bin_hits_exprs_down[,tall_exprs_name]),
#       na.rm=TRUE)
#y_min <- min(c(rpkm_left_up,gene_bin_hits_exprs_up[,tall_exprs_name],rpkm_left_down,gene_bin_hits_exprs_down[,tall_exprs_name]),
#       na.rm=TRUE)
y_min <- -6
y_max <- 6
par(mar=c(0.1,4,0.1,0.1))

# Add the unchanged in he middle of down and up
# Also, add the logFC of all genes as an additional control
boxplot(logFC_down, logFC_unchanged, logFC_all_genes,beside=TRUE,  notch=TRUE,col=c("blue","grey","dimgrey"),las=2,
        ylab="RNA logFC",na.action="na.omit", ylim=c(y_min,y_max),byt="n",pch=20,font=2,lwd=8,outline=FALSE)
#boxplot(logFC_down, logFC_up, beside=TRUE, notch=TRUE,col=c("blue","darkred"),na.action="na.omit", add=TRUE,las=2, ylim=c(y_min,y_max),bty="n",pch=20,lwd=8)
abline(h=0, lwd=10,col="black")

stripchart(list(logFC_down=logFC_down,logFC_unchanged=logFC_unchanged, logFC_all_genes=logFC_all_genes),
        col=c("blue","grey","dimgrey"),add=TRUE,pch=20,cex=1.4,vertical=TRUE,method="jitter")

boxplot(logFC_down, logFC_unchanged, logFC_all_genes,beside=TRUE, notch=TRUE,col=c("blue","grey","dimgrey"),na.action="na.omit", add=TRUE,las=2, ylim=c(y_min,y_max),bty="n",pch=20,lwd=8,outline=FALSE)
box(lwd=8)
axis(side = 2, lwd = 8,labels=FALSE)
#legend("topright", legend=c("Less TAD activity","Higher TAD activity"), fill=c("blue","darkred"))
#boxplot(rpkm_left_down,gene_bin_hits_exprs_down[,tall_exprs_name],beside=TRUE, main="Less TAD activity in T-ALL",notch=TRUE,col="lightgray",las=2,
#       names=box_names,ylab="logFPKM",ylim=c(y_min,y_max))
dev.off()

}


