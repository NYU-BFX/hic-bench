#!/usr/bin/Rscript
#$ -S /usr/bin/Rscript




################################
##### MAIN   ###################
################################


args <- commandArgs(trailingOnly=T)
if (length(args)!=3) {
  cat("USAGE: update-filtered-stats.r STATS-TSV N-INTRA-UNIQ N-INTER-UNIQ\n")
  quit(save="no")
}

stats = read.table(args[1],header=F,row.names=1,check.names=F)
n_intra_uniq = as.double(args[2])
n_inter_uniq = as.double(args[3])

stats["ds-duplicate-intra",1] = stats["ds-accepted-intra",1]-n_intra_uniq
stats["ds-duplicate-inter",1] = stats["ds-accepted-inter",1]-n_inter_uniq

stats["ds-accepted-intra",1] = n_intra_uniq
stats["ds-accepted-inter",1] = n_inter_uniq

stats[,2] = round(stats[,1]/stats["read-pairs",1],4)

options(scipen = 999)
write.table(cbind(rownames(stats),round(stats[,1],0),round(stats[,2],4)),sep='\t',quote=F,col.names=F,row.names=F)



