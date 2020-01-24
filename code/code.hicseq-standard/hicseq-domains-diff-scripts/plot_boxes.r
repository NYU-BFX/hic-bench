#R

args <- commandArgs()

less_activity_file <- args[3]
more_activity_file <- args[4]
unchanged_activity_file <- args[5]
out <- args[6]

less_activity <- read.csv(file=less_activity_file, header=FALSE,sep="\t")
more_activity <- read.csv(file=more_activity_file, header=FALSE,sep="\t")
unchanged_activity <- read.csv(file=unchanged_activity_file, header=FALSE,sep="\t")

results <- t.test(less_activity$V1, more_activity$V1)
results_less_unchanged <- t.test(less_activity$V1, unchanged_activity$V1, alternative = "less")
results_more_unchanged <- t.test(more_activity$V1, unchanged_activity$V1, alternative = "greater")
results_less <- t.test(less_activity$V1, alternative = "less")
results_more <- t.test(more_activity$V1, alternative = "greater")
results_unchanged <- t.test(unchanged_activity$V1)

matrix_results <- rbind(c(abs(results$estimate[1] - results$estimate[2]),results$p.value),
	c(abs(results_less_unchanged$estimate[1] - results_less_unchanged$estimate[2]),results_less_unchanged$p.value),
	c(abs(results_more_unchanged$estimate[1] - results_more_unchanged$estimate[2]),results_more_unchanged$p.value),
	c(results_less$estimate[1],results_less$p.value),
	c(results_more$estimate[1],results_more$p.value),
	c(results_unchanged$estimate[1],results_unchanged$p.value))
colnames(matrix_results) <- c("Mean diff","p-value")
rownames(matrix_results) <- c("less vs. more","less vs. unchanged", "more vs. unchanged", "less TAD act vs. 0", "more TAD act vs. 0", "unchanged TAD act vs. 0")

write.table(file=paste(out,"_ttest_results.tsv",sep=""),matrix_results,quote=FALSE,sep="\t")

mean_diff <- abs(results$estimate[1] - results$estimate[2])

#png(file=paste(out,"_boxplot.png",sep=""), width=2048, height=2048,pointsize=40)
#boxplot(less_activity$V1, more_activity$V1, col="white",outline=FALSE,las=2,lwd=4,names=c("less activity","more activity"),ylab="Aggregated CTCF fold-changes",
#	main=paste("Mean diff=",mean_diff,". p-val=",results$p.value,sep=""))
#stripchart(list(less_activity=less_activity$V1,more_activity=more_activity$V1),
#	col=c("blue","red"),add=TRUE,pch=20,cex=1.4,vertical=TRUE,method="jitter")
#abline(h=0,col="black",lwd=5)
#dev.off()


png(file=paste(out,"_boxplot.png",sep=""),width=1200,height=2048,pointsize=110)
#box_names <- c(case_right,case_left)
#y_max <- max(c(rpkm_left_up,gene_bin_hits_exprs_up[,tall_exprs_name],rpkm_left_down,gene_bin_hits_exprs_down[,tall_exprs_name]),
#       na.rm=TRUE)
#y_min <- min(c(rpkm_left_up,gene_bin_hits_exprs_up[,tall_exprs_name],rpkm_left_down,gene_bin_hits_exprs_down[,tall_exprs_name]),
#       na.rm=TRUE)
#y_min <- -16
#y_max <- 20
par(mar=c(0.1,4,0.1,0.1))
boxplot(less_activity$V1,unchanged_activity$V1, more_activity$V1, beside=TRUE,  notch=TRUE,col=c("blue","grey","darkred"),las=2,
        ylab="Aggregated ATAC seq fold-changes",na.action="na.omit",byt="n",pch=20,font=2,lwd=8,axis=FALSE,outline=FALSE)
abline(h=0, lwd=10,col="black")
stripchart(list(less_activity=less_activity$V1,unchanged_activity=unchanged_activity$V1,more_activity=more_activity$V1),
	col=c("blue","grey","darkred"),add=TRUE,pch=20,cex=1.4,vertical=TRUE,method="jitter")
boxplot(less_activity$V1, unchanged_activity$V1, more_activity$V1, beside=TRUE, notch=TRUE,col=c("blue","grey","darkred"),na.action="na.omit", add=TRUE,las=2,bty="n",pch=20,lwd=8,axis=FALSE,
	outline=FALSE)
box(lwd=8)
axis(side = 2, lwd = 8,labels=FALSE)
#legend("topright", legend=c("Less TAD activity","Higher TAD activity"), fill=c("blue","darkred"))
#boxplot(rpkm_left_down,gene_bin_hits_exprs_down[,tall_exprs_name],beside=TRUE, main="Less TAD activity in T-ALL",notch=TRUE,col="lightgray",las=2,
#       names=box_names,ylab="logFPKM",ylim=c(y_min,y_max))
dev.off()



stop(0)

results <- cor.test(less_activity$V1, less_activity$V2,method="spearman")
png(file=paste(out,"_scatter-Tcell-boundary.png",sep=""),width=2048, height=2048,pointsize=40)
plot(less_activity$V1, less_activity$V2,pch=20,col="black",cex=1.5,xlab="CTCF Fold",ylab="TAD boundary - intra-TAD mean diff",main=paste("R=",results$estimate,sep=""))
dev.off()

results <- cor.test(more_activity$V1, more_activity$V2,method="spearman")
png(file=paste(out,"_scatter-more_activity-boundary.png",sep=""),width=2048, height=2048,pointsize=40)
plot(more_activity$V1, more_activity$V2,pch=20,col="black",cex=1.5,xlab="CTCF Fold",ylab="TAD boundary - intra-TAD mean diff",main=paste("R=",results$estimate,sep=""))
dev.off()
