### Functions ##############################################################################
getEigenMat=function(pca1.matrix,centrotelo){
  df.sizes=data.frame(sampleName=names(pca1.matrix),compA.size=NA,compB.size=NA,compA.meanPC1=NA,compB.meanPC1=NA,compA.medianPC1=NA,compB.medianPC1=NA,stringsAsFactors = F) #to store compartments size info
  centrotelo.gr=makeGRangesFromDataFrame(centrotelo)
  rownames(pca1.matrix)=gsub(rownames(pca1.matrix),pattern = "-",replacement = ":")
  t=unlist(strsplit(rownames(pca1.matrix),split = ":"))
  pca1.matrix$chr=t[seq(1,length(t),3)]
  pca1.matrix$start=t[seq(2,length(t),3)]
  pca1.matrix$end=t[seq(3,length(t),3)]
  pca1.matrix.gr=makeGRangesFromDataFrame(pca1.matrix,keep.extra.columns = T)
  pca1.matrix=as.data.frame(subsetByOverlaps(pca1.matrix.gr,centrotelo.gr,invert = T))
  names(pca1.matrix)[1]="chr"
  pca1.matrix=pca1.matrix[c(1:3,6:ncol(pca1.matrix))]
  pca1.matrix=pca1.matrix[order(pca1.matrix$chr,pca1.matrix$start),]
  pca1.matrix[,4:ncol(pca1.matrix)]=pca1.matrix[,4:ncol(pca1.matrix)]/100
  binsize=median(pca1.matrix$end-pca1.matrix$start,na.rm = T) # get bin size
  
  for (S1 in names(pca1.matrix)[4:ncol(pca1.matrix)]){
    sampleLabel=S1
    df.sizes$compA.size[df.sizes$sampleName==S1]=sum(pca1.matrix[,S1] > 0,na.rm = T)*binsize
    df.sizes$compB.size[df.sizes$sampleName==S1]=sum(pca1.matrix[,S1] < 0,na.rm = T)*binsize
    df.sizes$compA.meanPC1[df.sizes$sampleName==S1]=mean(pca1.matrix[pca1.matrix[,S1] > 0,S1],na.rm = T)
    df.sizes$compB.meanPC1[df.sizes$sampleName==S1]=mean(pca1.matrix[pca1.matrix[,S1] < 0,S1],na.rm = T)
    df.sizes$compA.medianPC1[df.sizes$sampleName==S1]=median(pca1.matrix[pca1.matrix[,S1] > 0,S1],na.rm = T)
    df.sizes$compB.medianPC1[df.sizes$sampleName==S1]=median(pca1.matrix[pca1.matrix[,S1] < 0,S1],na.rm = T)
    df.sizes$sampleName[df.sizes$sampleName==S1]=sampleLabel
  }
  df.sizes[,2:ncol(df.sizes)]=apply(df.sizes[,2:ncol(df.sizes)],2,as.numeric)
  df.sizes$perc.compA=df.sizes$compA.size/(df.sizes$compA.size+df.sizes$compB.size)
  df.sizes$perc.compB=df.sizes$compB.size/(df.sizes$compA.size+df.sizes$compB.size)
  df.sizes$size.diff.Mb=round(abs(df.sizes$compB.size-df.sizes$compA.size)/1000000,digits = 1)
  
  return(list("mat"=pca1.matrix,"df.sizes"=df.sizes))
}

getBarPlot=function(df.sizes, x, y, grid){
  df.sizes.gg=data.frame(size=c(df.sizes[,"compA.size"],
                                df.sizes[,"compB.size"]),
                         group=rep(df.sizes$sampleName,2),
                         class=rep(c("compartment A","compartment B"),each=nrow(df.sizes)),
                         perc.size=c(df.sizes$perc.compA,df.sizes$perc.compB))
  df.sizes.gg$sampleCompartment=paste0(df.sizes.gg$group,".",df.sizes.gg$class)
  df.sizes.gg$size=round(df.sizes.gg$size/1000000)
  if(grid==T){ df.sizes.gg=df.sizes.gg[as.character(df.sizes.gg$group)==names(mat)[x] | as.character(df.sizes.gg$group)==names(mat)[x],] }
  df.sizes.gg=df.sizes.gg[order(df.sizes.gg$group),]
  df.sizes.gg$sizeLab.pos=df.sizes.gg$size*0.5
  df.sizes.gg$perc.size=round(df.sizes.gg$perc.size*100,1)
  df.sizes.gg$sizeLab.pos[df.sizes.gg$class=="compartment A"]=(df.sizes.gg$size[df.sizes.gg$class=="compartment A"])/2+df.sizes.gg$size[df.sizes.gg$class=="compartment B"]
  df.sizes.gg$sizeLabel[df.sizes.gg$class=="compartment A"]=paste0(df.sizes.gg$size[df.sizes.gg$class=="compartment A"]," Mb\nA (",df.sizes.gg$perc.size[df.sizes.gg$class=="compartment A"],"%)")
  df.sizes.gg$sizeLabel[df.sizes.gg$class=="compartment B"]=paste0(df.sizes.gg$size[df.sizes.gg$class=="compartment B"]," Mb\nB (",df.sizes.gg$perc.size[df.sizes.gg$class=="compartment B"],"%)")
  
  ggplot(data = df.sizes.gg, aes(x = group, y = size, fill = class)) + 
    geom_bar(stat = "identity",colour="black",position = "stack",width = 0.6)+
    xlab("")+
    ylab("compartments size (Mb)")+
    scale_fill_manual(values=c("red","blue"))+
    geom_text(data = df.sizes.gg,
              aes(x = group,
                  y = sizeLab.pos,
                  label = sizeLabel),colour="white",
              size = 4,
              vjust = 0.5)+
    scale_x_discrete(labels= as.character(labels))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    guides(fill=guide_legend(""))
}

getDensPlot=function(mat, x, y){
  temp=data.frame(PC1=c(mat[,names(mat)[x]],mat[,names(mat)[y]]),group=rep(c(names(mat)[y],names(mat)[x]),each=nrow(mat)),stringsAsFactors = F)
  ggplot(data = temp, aes(x=PC1,fill=group))+
    geom_density(show.legend = T,alpha=0.5)+
    scale_fill_manual(values=c("green","purple"),labels=labels[c(x,y)])+
    theme(legend.position = c(.95, .95),
          legend.justification = c("right", "top"),
          legend.box.just = "right",
          legend.margin = margin(6, 6, 6, 6),
          legend.title = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+
    geom_vline(xintercept = 0,linetype = "dotdash",color="darkred")

}

getHexPlot <- function(mat, x, y, grid){
  if(grid==T){sz=3} else {sz=5}
  correl=round(cor(mat[,names(mat)[x]],mat[,names(mat)[y]]),digits = 3)
  ggplot(data = mat, aes_string(x=names(mat)[x],y=names(mat)[y]))+ 
    stat_binhex(bins=70) +
    xlab(labels[x])+
    ylab(labels[y])+
    scale_fill_gradientn(colours=c("blue","orange","red"),trans="log10","  log10 count\n(eigenvector-1)\n")+
    geom_abline(slope = 1,intercept = 0,linetype=2)+
    annotate(geom = 'text', label = paste0('r = ',correl), x = -Inf, y = Inf, hjust = -0.2, vjust = 1.5,size=sz)+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))
}

getEigenPairs=function(mat){
  grid=T
  p=ggpairs(mat,lower="blank",diag="blank",upper="blank",columnLabels = labels,legend=c(2,1))+
    theme_bw(base_size=14)+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.direction = "vertical",
          legend.position = "none")+
    theme(axis.text=element_text(size=8,color="black"),
          axis.title=element_text(size=8,face="bold"),
          strip.background = element_rect(fill="grey"),
          legend.position = "left")
  seq=1:ncol(mat)
  for (x in seq){
    for (y in seq){
      if (y>x){
        p=putPlot(p, getHexPlot(mat=mat,x=x,y=y,grid=T),y,x)
        p=putPlot(p, getDensPlot(mat=mat,x=x,y=y),x,y)
      }
    }
  }
  for (x in seq){ p=putPlot(p, getBarPlot(df=df.sizes,x=x,y=x,grid=T),x,x) }
  return(print(p))
}

getEigenPCA=function(mat,df.color){
  pca = prcomp(t(mat))
  names = colnames(mat)
  #fac = factor(sapply(names,function(x){strsplit(x,'.')[[1]][1]}))
  #colours = rep(c(brewer.pal(7,"Set1"),brewer.pal(7,"Set2"),brewer.pal(7,"Set3")),nlevels(fac))[1:nlevels(fac)]
  #groupColor=data.frame(group=as.character(unique(fac)),color=as.character(colours))
  #df=data.frame(sampleName=names,group=fac,color=NA)
  #for (group in df$group){df$color[df$group==group]=as.character(groupColor$color[groupColor$group==group])}
  print(autoplot(pca,label = F,label.size = 2,size=5,col=df$color)+
          geom_text_repel(aes(label=labels),
                          segment.size = 0.5,
                          force = 3,size = 4,
                          point.padding = 1,
                          segment.alpha = 1)+
          geom_point(size=2)+
          theme(panel.background = element_rect(fill = "white"),
                panel.border = element_rect(colour = "black",fill=NA,size = 2),
                axis.text.x = element_text(size = 20,colour = "black"),
                axis.text.y = element_text(size = 20,colour = "black"),  
                axis.title.x = element_text(size = 20,colour = "black"),
                axis.title.y = element_text(size = 20,colour = "black"),
                legend.position="none"))
}

getEigenHmap=function(mat,tot.bins=3000){ #tot.bins = number of top most variable pc1 bins to plot
  n_samples=ncol(mat)
  mat$sd=rowSds(as.matrix(mat),na.rm=T)
  mat=mat[order(mat$sd,decreasing = T),]
  mat=mat[,1:n_samples]
  names=names(mat)
  fac = factor(sapply(names,function(x){strsplit(x,'.')[[1]][1]}))
 
  heatmap.2(as.matrix(mat[1:tot.bins,1:n_samples]),
            col=colorRampPalette(rev(brewer.pal(n = 10, name ="RdBu")))(100),
            dendrogram = "both",
            scale= "row",
            labCol=labels,
            hclustfun = function(x){hclust(x, method="ward.D2")},
            distfun = function(x) {dist(x, method="manhattan")},
            labRow = FALSE,
            cexRow=1,
            cexCol=1,
            margins=c(12,7),
            trace="none",
            srtCol=45,
            key = TRUE,
            keysize = 1.5,
            density.info=c("none"))
}

# classify switched pc1 bins (pairwise-comparison)
classifyPc1Bins=function(mat){
  df.switch=data.frame(sampleName1=character(choose(ncol(mat),2)),sampleName2=NA,AA.number=NA,AB.number=NA,BB.number=NA,BA.number=NA,stringsAsFactors = F)
  seq=1:ncol(mat)
  nr=1
  for(x in seq){
    for (y in seq){
      if(y>x){
        df.mat=mat[,c(names(mat)[x],names(mat)[y])]
        df.mat$pc1.diff=df.mat[,1]-df.mat[,2]
        pc1Diff_q50=quantile(abs(df.mat$pc1.diff),0.5)
        #pc1.diff.cut=pc1Diff_q50
        pc1.diff.cut=0.0008
        df.mat$switch=NA
        df.mat$switch[df.mat[,1] > 0 & df.mat[,2,] > 0]="AA"
        df.mat$switch[df.mat[,1] < 0 & df.mat[,2,] < 0]="BB"
        df.mat$switch[df.mat[,1] > 0 & df.mat[,2,] < 0 & abs(df.mat$pc1.diff) >= pc1.diff.cut]="AB"
        df.mat$switch[df.mat[,1] < 0 & df.mat[,2,] > 0 & abs(df.mat$pc1.diff) >= pc1.diff.cut]="BA"
        df.mat$switch[df.mat[,1] > 0 & df.mat[,2,] < 0 & abs(df.mat$pc1.diff) < pc1.diff.cut]="AA"
        df.mat$switch[df.mat[,1] < 0 & df.mat[,2,] > 0 & abs(df.mat$pc1.diff) < pc1.diff.cut]="BB"
        df.switch$AA.number[nr]=sum(df.mat$switch=="AA",na.rm = T)
        df.switch$AB.number[nr]=sum(df.mat$switch=="AB",na.rm = T)
        df.switch$BB.number[nr]=sum(df.mat$switch=="BB",na.rm = T)
        df.switch$BA.number[nr]=sum(df.mat$switch=="BA",na.rm = T)
        df.switch$sampleName1[nr]=labels[y]
        df.switch$sampleName2[nr]=labels[x]
        write.table(df.mat,paste0(out_dir,'/pc1.switch_',names(mat)[y],"_vs_",names(mat)[x],'.tsv'),sep='\t',row.names = F,col.names=T,quote=F)
        nr=nr+1
      }
    }
  }
  df.switch$comparisonID=paste0(df.switch$sampleName1,":",df.switch$sampleName2)
  return(df.switch)
}

#barplots switched-pc1-bins
getBarSwitch=function(df.switch){
  t=as.data.frame(t(df.switch[,c("AA.number","AB.number","BB.number","BA.number")]))
  names(t)=df.switch$comparisonID
  df.t=data.frame(number=numeric(4*ncol(t)),class=NA,comparisonID=NA,stringsAsFactors = F)
  
  r=1
  for (i in 1:ncol(t)){
    df.t$number[r:(r+3)]=t[,i]
    df.t$comparisonID[r:(r+3)]=names(t)[i]
    df.t$class[r:(r+3)]=row.names(t)
    r=r+4
  }
  
  df.t$class=gsub(df.t$class,pattern = ".number",replacement = "")
  df.t$comparisonID=gsub(df.t$comparisonID,pattern = ":",replacement = " vs ",fixed = T)
  
  print(ggplot(data = df.t, aes(x = comparisonID, y = number, fill = class)) + 
          geom_bar(stat = "identity",colour="black",position = "stack",width = 0.6,size=0.5)+
          xlab("")+
          ylab("number of eigenvector-1 bins")+
          scale_fill_manual(values=c("darkred","blue","red","darkblue"))+
          coord_flip()+
          guides(fill=guide_legend("")))
}

### RUN ################################################################################

# load libraries
library("hexbin")
library("lattice")
library("RColorBrewer")
library("ggplot2")
library("GGally")
library("matrixStats")
library("gplots")
library("ggfortify")
library("ggrepel")
library("pheatmap")
library("ggforce")
library("GenomicRanges")

# process command-line arguments
cmdline_args <- commandArgs(trailingOnly=T);

for (p in c("optparse")) # install packages
  if (!suppressPackageStartupMessages(require(p,character.only=TRUE,quietly=TRUE,warn.conflicts=FALSE))) {
    install.packages(p,repos="http://cran.rstudio.com/") 
    library(p,character.only=TRUE,quietly=TRUE,verbose=FALSE)
  }

option_list <- list(
  make_option(c("-v","--verbose"), action="store_true",default=FALSE, help="Print more messages."),
  make_option(c("-o","--output-dir"), default="", help="Output directory (required) [default \"%default\"]."),
  make_option(c("-L","--sample-labels"), default="", help="A file containing sample labels in its first column [default \"%default\"]."),
  make_option(c("-c","--centrotelo-file"), default="", help="A bed file containing centromere and telomere coordinates [default \"%default\"]."),
  make_option(c("-m","--pca1-matrix"), default="", help="A matrix containing the samples PC1 values for each sample [default \"%default\"]."),
  make_option(c("--show-text"), action="store_true",default=TRUE, help="Display sample label text on PCA plot."),
  make_option(c("--use-short-names"), action="store_true",default=FALSE, help="Use only second part of the name (after the colon)."),
  make_option(c("--plain"), action="store_true",default=FALSE, help="Create only PCA and heatmap plots.")
);
usage = 'scripts-compartments-stats.r [OPTIONS] MATRIX';

# get command line options & input arguments
arguments <- parse_args(args=cmdline_args, OptionParser(usage=usage,option_list=option_list), positional_arguments=c(0,Inf));
opt <- arguments$options;
files <- opt$'pca1-matrix'
if (length(files)!=1) { write('Error: this operation requires an input matrix!',stderr()); quit(save='no'); }
out_dir = opt$'output-dir'
sample_labels = opt$'sample-labels'
show_text = opt$'show-text'
use_short_names = opt$'use-short-names'
plain = opt$'plain'
centrotelo_file = opt$'centrotelo-file'

# create output directories
if (file.exists(out_dir)==FALSE) { dir.create(out_dir) } else { write('Warning: output directory already exists, it will be overwritten!',stderr()) }
if (file.exists(paste0(out_dir,"/hexbin_plots"))==FALSE) { dir.create(paste0(out_dir,"/hexbin_plots")) } else { write('Warning: output directory already exists, it will be overwritten!',stderr()) }
if (file.exists(paste0(out_dir,"/density_plots"))==FALSE) { dir.create(paste0(out_dir,"/density_plots")) } else { write('Warning: output directory already exists, it will be overwritten!',stderr()) }

# load input files
pca1.matrix=read.table(files,stringsAsFactors = F) 
centrotelo=read.table(centrotelo_file,stringsAsFactors = F) 
bed.cols=c("chr","start","end")
names(centrotelo)=bed.cols

matAndSizes=getEigenMat(pca1.matrix=pca1.matrix,centrotelo=centrotelo) #get pca1 matrix and compartment sizes
mat=as.data.frame(matAndSizes$mat)
df.sizes=as.data.frame(matAndSizes$df.sizes)
mat=mat[,4:ncol(mat)]

# set label format
if (sample_labels!="") { colnames(mat) = t(read.table(sample_labels,header=F)[,1]) } 
labels=names(mat) 
test.1st.char=as.numeric(substr(x = labels,start = 1,stop = 1)) #to avoid ggplot error when first character of label is numeric
test.1st.char=test.1st.char[!is.na(test.1st.char)]
if(length(test.1st.char)>0){labels=paste0("S",labels)}

# set color groups
fac = factor(sapply(labels,function(x){strsplit(x,':',fixed = T)[[1]][1]}))
colours = rep(c(brewer.pal(7,"Set1"),brewer.pal(7,"Set2"),brewer.pal(7,"Set3")),nlevels(fac))[1:nlevels(fac)]
groupColor=data.frame(group=as.character(unique(fac)),color=as.character(colours))
df=data.frame(sampleName=labels,group=fac,color=NA)
for (group in df$group){df$color[df$group==group]=as.character(groupColor$color[groupColor$group==group])}

# avoid ggplot errors when using certain symbols
labels=gsub(labels,pattern = "-",replacement = "_")
labels=gsub(labels,pattern = ":",replacement = ".")
colnames(mat)=labels
df.sizes$sampleName=labels #for labels consistency
short_names = as.vector(sapply(labels,function(x){strsplit(x,'.',fixed=T)[[1]][2]}))
if (show_text) { if (use_short_names) { labels = short_names } else { labels = labels } } else { labels = NULL }


### MAKE PLOTS ### 
if(length(names(mat))>2){ # make pca & hmap only when the number of samples > 2
  pdf(paste(out_dir,'/eigenPCA.pdf',sep=''),useDingbats=FALSE)
  getEigenPCA(mat,df.color = df.color)
  dev.off()
  
  pdf(paste(out_dir,'/eigenHMAP.pdf',sep=''),width =10,useDingbats=FALSE)
  getEigenHmap(mat)
  dev.off()
}

if (plain==FALSE){
  # make hexbin-bar-density grid plot (all vs all)
  pdf(paste(out_dir,'/eigenPairsGrid.pdf',sep=''),height = 10,width =14,useDingbats=FALSE)
  getEigenPairs(mat)
  dev.off()
  
  # get separate hexbinplots (all vs all)
  seq=1:ncol(mat)
  for(x in seq){
    for (y in seq){
      if(y>x){
        pdf(paste(out_dir,paste0('/hexbin_plots/eigenHexbinPlot_',names(mat)[y],"_vs_",names(mat)[x],'.pdf'),sep=''),useDingbats=FALSE)
        print(getHexPlot(mat,x,y,grid=F))
        dev.off()
        
        pdf(paste(out_dir,paste0('/density_plots/eigenDensityPlot_',names(mat)[y],"_vs_",names(mat)[x],'.pdf'),sep=''),useDingbats=FALSE)
        print(getDensPlot(mat,x,y))
        dev.off()
      }
    }
  }
  
  # get all barplots in one plot
  pdf(paste0(out_dir,'/size_barplots.pdf'),width = 9,useDingbats=FALSE)
  print(getBarPlot(df.sizes,x=1,y=1,grid = F))
  dev.off()
  
  # classify pc1 bin in AA, AB, BB, BA
  df.switch=classifyPc1Bins(mat)

    #plot bars for pc1 classified bins
  pdf(paste0(out_dir,'/switch_bins_barplot.pdf'),width = 12,useDingbats=FALSE)
  getBarSwitch(df.switch)
  dev.off()
}
