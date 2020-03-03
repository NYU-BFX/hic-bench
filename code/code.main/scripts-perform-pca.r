#!/usr/bin/env Rscript
#$ -S /usr/bin/env Rscript

##
## USAGE: scripts-perform-pca.r [OPTIONS] MATRIX
##
plotPCA <- function(mat,show_text,use_short_names,plain){
  
  pca = prcomp(t(mat))
  
  if (plain==FALSE) {heatmap.2(pca$x,scale='column',key=FALSE,trace='none',margins=c(1,5),Colv=FALSE,dendrogram='row',labCol=NULL,cexRow=0.5)}
  
  names = colnames(mat)
  fac = factor(sapply(names,function(x){strsplit(x,':')[[1]][1]}))
  short_names = as.vector(sapply(names,function(x){strsplit(x,':')[[1]][2]}))
  if (show_text) { if (use_short_names) { labels = short_names } else { labels = names } } else { labels = NULL }
  colours = rep(c(brewer.pal(7,"Set1"),brewer.pal(7,"Set2"),brewer.pal(7,"Set3")),nlevels(fac))[1:nlevels(fac)]
  
  if (plain==TRUE) { ## Modified by Javier ##
  print(autoplot(pca,label = F,label.size = 2,colour=colours,size=5)+
          geom_text_repel(aes(label=colnames(mat)),
                          segment.size = 0.2,
                          force = 3,size=5,
                          point.padding=2,
                          segment.alpha=0)+
          scale_color_manual(values=colours)+
          geom_point(size=2)+
          theme(panel.background = element_rect(fill = "white"),
                panel.border = element_rect(colour = "black",fill=NA,size = 2),
                axis.text.x = element_text(size = 20,colour = "black"),
                axis.text.y = element_text(size = 20,colour = "black"),  
                axis.title.x = element_text(size = 20,colour = "black"),
                axis.title.y = element_text(size = 20,colour = "black")))
    }
  
  return(pca)
}




################################
##### MAIN   ###################
################################


cmdline_args <- commandArgs(trailingOnly=T);

# install packages
for (p in c("optparse")) 
  if (!suppressPackageStartupMessages(require(p,character.only=TRUE,quietly=TRUE,warn.conflicts=FALSE))) {
    install.packages(p,repos="http://cran.rstudio.com/") 
    library(p,character.only=TRUE,quietly=TRUE,verbose=FALSE)
  }

# process command-line arguments
option_list <- list(
  make_option(c("-v","--verbose"), action="store_true",default=FALSE, help="Print more messages."),
  make_option(c("-o","--output-dir"), default="", help="Output directory (required) [default \"%default\"]."),
  make_option(c("-L","--sample-labels"), default="", help="A file containing sample labels in its first column [default \"%default\"]."),
  make_option(c("--show-text"), action="store_true",default=FALSE, help="Display sample label text on PCA plot."),
  make_option(c("--use-short-names"), action="store_true",default=FALSE, help="Use only second part of the name (after the colon)."),
  make_option(c("--plain"), action="store_true",default=FALSE, help="Create only one page (PC1 vs PC2).")
);
usage = 'perform_pca.r [OPTIONS] MATRIX';

# get command line options (if help option encountered print help and exit)
arguments <- parse_args(args=cmdline_args, OptionParser(usage=usage,option_list=option_list), positional_arguments=c(0,Inf));
opt <- arguments$options;
files <- arguments$args;
if (length(files)!=1) { write('Error: this operation requires an input matrix!',stderr()); quit(save='no'); }

# input arguments
out_dir = opt$'output-dir'
sample_labels = opt$'sample-labels'
show_text = opt$'show-text'

# create output directory
if (out_dir=="") { write('Error: please specify output directory!',stderr()); quit(save='no'); }
if (file.exists(out_dir)==FALSE) { dir.create(out_dir) } else { write('Warning: output directory already exists, it will be overwritten!',stderr()) }

# load data
if (opt$'verbose'==TRUE) write('Loading input matrix...',stderr())
x = read.table(files[1],stringsAsFactors = F)
if (sample_labels!="") { colnames(x) = t(read.table(sample_labels,header=F)[,1]) }


# load libraries
library("RColorBrewer")
library("lattice")
library("gplots")
library(ggfortify)
library(ggrepel)

# PCA on raw input matrix
if (opt$'verbose'==TRUE) write('Performing PCA on input matrix...',stderr())
pdf(paste(out_dir,'/report.raw.pdf',sep=''),useDingbats=FALSE)
x_pca = plotPCA(x,show_text=show_text,use_short_names=opt$"use-short-names",plain=opt$"plain")
dev.off()

# PCA on mean-normalized matrix
if (opt$'verbose'==TRUE) write('Performing PCA on mean-normalized matrix...',stderr())
z = t(t(x)/apply(x,2,mean))
pdf(paste(out_dir,'/report.mnorm.pdf',sep=''),useDingbats=FALSE)
z_pca = plotPCA(z,show_text=show_text,use_short_names=opt$"use-short-names",plain=opt$"plain")
dev.off()

# PCA on quantile-normalized matrix
if (opt$'verbose'==TRUE) write('Performing PCA on quantile-normalized matrix...',stderr())
suppressMessages(library("preprocessCore"))
y = normalize.quantiles(as.matrix(x))
rownames(y) = rownames(x)
colnames(y) = colnames(x)
pdf(paste(out_dir,'/report.qnorm.pdf',sep=''),useDingbats=FALSE)
y_pca = plotPCA(y,show_text=show_text,use_short_names=opt$"use-short-names",plain=opt$"plain")
dev.off()

#if (opt$'verbose'==TRUE) write('Saving results...',stderr())
#save(x_pca,y_pca,file=paste(out_dir,'/data.RData',sep=''))

if (opt$'verbose'==TRUE) write('Done.',stderr())
quit(save='no')




