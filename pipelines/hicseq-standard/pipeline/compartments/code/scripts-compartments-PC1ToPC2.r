### USAGE ###
# EXAMPLE:
# Rscript ./code/scripts-compartments-PC1ToPC2.r ./results/compartments.by_group.homer.res_100kb/filter.by_sample.mapq_20/align.by_sample.bowtie2/fadu-WT/ chr3,chr17

argv = commandArgs(trailingOnly = TRUE)
samplePath = argv[1L] #path to the sample/group results directory
chromsList = argv[2L] #comma separeted list of the chromosomes to be flipped

chromsToFlip=unlist(strsplit(chromsList,split = ","))
head(chromsToFlip)
bedgraph=read.csv(paste0(samplePath,"/pca_HKgenesFix.PC1.bedGraph"),sep="\t",header = F,skip=1,,stringsAsFactors = F,col.names = c("chr","start","end","pc1"))
bedgraph$bin.ID=paste(bedgraph$chr,bedgraph$start,bedgraph$end,sep = ":")
head(bedgraph)
bedgraph.PC2=read.csv(paste0(samplePath,"/pca_HKgenesFix.PC1.PC2.txt"),sep="\t",header = T,stringsAsFactors = F)
bedgraph.PC2$bin.ID=paste(bedgraph.PC2$chr,bedgraph.PC2$start,bedgraph.PC2$end,sep = ":")
head(bedgraph.PC2)
bedgraph.PC2=bedgraph.PC2[bedgraph.PC2$bin.ID %in% bedgraph$bin.ID,]
print("after")
write.table(bedgraph,paste0(samplePath,"/pca_HKgenesFix.PC1_old.bedGraph"),quote = F,col.names = F,row.names = F,sep="\t")

# Flip to PC2 values in the selected chromosomes
bedgraph$pc1[bedgraph$chr %in% chromsToFlip]=bedgraph.PC2$PC2[bedgraph.PC2$bin.ID %in% bedgraph$bin.ID[bedgraph$chr %in% chromsToFlip]]
write.table(bedgraph,paste0(samplePath,"/pca_HKgenesFix.PC1.bedGraph"),quote = F,col.names = F,row.names = F,sep="\t")

print("Done!")
print(paste0("The PC1 values of ",chromsList," has been replaced by the PC2 values"))
print("The 'pca_HKgenesFix.PC1.bedGraph' file has been replaced with a version that incorporates those changes")
print("The previous 'pca_HKgenesFix.PC1.bedGraph' file has been saved as 'pca_HKgenesFix.PC1_old.bedGraph'")
print("The 'pca_HKgenesFix_Acompartments.bed' file has been replaced with a version that incorporates those changes")
print("If you already ran the compartments-stats step, please re-run it to get the updated stats")
