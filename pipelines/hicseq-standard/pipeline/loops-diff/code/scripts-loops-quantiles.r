argv = commandArgs(trailingOnly = TRUE)
infile = argv[1L]
outdir = argv[2L]
print(infile)
library(ggplot2)

label="loops_filtered_cpm"
options(scipen=10000)
x=read.table(infile,header = F,stringsAsFactors = F)
names(x)[7]="cpm"
x$distance=x$V5-x$V3

qdist=quantile(x$distance,c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),na.rm=T)
qdist1=x[x$distance<=qdist[1],]
qdist2=x[x$distance>qdist[1] & x$cpm<=qdist[2],]
qdist3=x[x$distance>qdist[2] & x$cpm<=qdist[3],]
qdist4=x[x$distance>qdist[3] & x$cpm<=qdist[4],] 
qdist5=x[x$distance>qdist[4] & x$cpm<=qdist[5],] 
qdist6=x[x$distance>qdist[5] & x$cpm<=qdist[6],] 
qdist7=x[x$distance>qdist[6] & x$cpm<=qdist[7],] 
qdist8=x[x$distance>qdist[7] & x$cpm<=qdist[8],] 
qdist9=x[x$distance>qdist[8] & x$cpm<=qdist[9],] 
qdist10=x[x$distance>qdist[9],] 
  
write.table(qdist1,paste0(outdir,label,"_qdist1.bedpe"),sep = "\t",quote=F,row.names = F,col.names = F)
write.table(qdist2,paste0(outdir,label,"_qdist2.bedpe"),sep = "\t",quote=F,row.names = F,col.names = F)
write.table(qdist3,paste0(outdir,label,"_qdist3.bedpe"),sep = "\t",quote=F,row.names = F,col.names = F)
write.table(qdist4,paste0(outdir,label,"_qdist4.bedpe"),sep = "\t",quote=F,row.names = F,col.names = F)
write.table(qdist5,paste0(outdir,label,"_qdist5.bedpe"),sep = "\t",quote=F,row.names = F,col.names = F)
write.table(qdist6,paste0(outdir,label,"_qdist6.bedpe"),sep = "\t",quote=F,row.names = F,col.names = F)
write.table(qdist7,paste0(outdir,label,"_qdist7.bedpe"),sep = "\t",quote=F,row.names = F,col.names = F)
write.table(qdist8,paste0(outdir,label,"_qdist8.bedpe"),sep = "\t",quote=F,row.names = F,col.names = F)
write.table(qdist9,paste0(outdir,label,"_qdist9.bedpe"),sep = "\t",quote=F,row.names = F,col.names = F)
write.table(qdist10,paste0(outdir,label,"_qdist10.bedpe"),sep = "\t",quote=F,row.names = F,col.names = F)
  
qdist.df=as.data.frame(qdist)
write.csv(qdist.df,paste0(outdir,label,"_qdist.csv"),row.names = F)
  
qcpm10=quantile(x$cpm,c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),na.rm=T)
q1=x[x$cpm<=qcpm10[1],] 
q2=x[x$cpm>qcpm10[1] & x$cpm<=qcpm10[2],]
q3=x[x$cpm>qcpm10[2] & x$cpm<=qcpm10[3],]
q4=x[x$cpm>qcpm10[3] & x$cpm<=qcpm10[4],]
q5=x[x$cpm>qcpm10[4] & x$cpm<=qcpm10[5],]
q6=x[x$cpm>qcpm10[5] & x$cpm<=qcpm10[6],]
q7=x[x$cpm>qcpm10[6] & x$cpm<=qcpm10[7],]
q8=x[x$cpm>qcpm10[7] & x$cpm<=qcpm10[8],]
q9=x[x$cpm>qcpm10[8] & x$cpm<=qcpm10[9],] 
q10=x[x$cpm>qcpm10[9],]
  
write.table(q1,paste0(outdir,label,"_q1.bedpe"),sep = "\t",quote=F,row.names = F,col.names = F)
write.table(q2,paste0(outdir,label,"_q2.bedpe"),sep = "\t",quote=F,row.names = F,col.names = F)
write.table(q3,paste0(outdir,label,"_q3.bedpe"),sep = "\t",quote=F,row.names = F,col.names = F)
write.table(q4,paste0(outdir,label,"_q4.bedpe"),sep = "\t",quote=F,row.names = F,col.names = F)
write.table(q5,paste0(outdir,label,"_q5.bedpe"),sep = "\t",quote=F,row.names = F,col.names = F)
write.table(q6,paste0(outdir,label,"_q6.bedpe"),sep = "\t",quote=F,row.names = F,col.names = F)
write.table(q7,paste0(outdir,label,"_q7.bedpe"),sep = "\t",quote=F,row.names = F,col.names = F)
write.table(q8,paste0(outdir,label,"_q8.bedpe"),sep = "\t",quote=F,row.names = F,col.names = F)
write.table(q9,paste0(outdir,label,"_q9.bedpe"),sep = "\t",quote=F,row.names = F,col.names = F)
write.table(q10,paste0(outdir,label,"_q10.bedpe"),sep = "\t",quote=F,row.names = F,col.names = F)
qcpm10.df=as.data.frame(qcpm10)
write.csv(qcpm10.df,paste0(outdir,label,"_qcpm.csv"),row.names = F)
