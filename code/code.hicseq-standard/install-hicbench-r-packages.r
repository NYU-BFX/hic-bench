#!/usr/bin/env Rscript

write("Checking installation of R packages", stderr())


packages = c("igraph", "flsa", "genlasso", "optparse", "ggplot2", "pastecs", "plotrix",
             "reshape2", "zoo", "plyr", "data.table", "gridExtra", "scales", "grid",
             "RColorBrewer", "stringr", "chron", "corrplot", 'dplyr', 'tidyr', 'stringr', 'lubridate', 'readr', 'cowplot', 'pheatmap', 'tools', 
             "ggrepel", "ggfortify", "argparser")

for (p in packages)
  if (!suppressPackageStartupMessages(require(p, character.only=TRUE, quietly=TRUE, warn.conflicts=FALSE))) {
    install.packages(p, repos="http://cran.us.r-project.org")
  }
  
#install.packages('code/hicrep_1.0.1.tar.gz')

