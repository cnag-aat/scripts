#!/usr/bin/env Rscript
require(treemap)
library(RColorBrewer)
argv <- commandArgs(trailingOnly = TRUE)
blucol<-brewer.pal(9, "Blues")
filename<-argv[1]
plottitle<-argv[2]
assembly<-read.table(filename)
pdf(paste0(filename,".treemap.pdf"))
treemap(assembly,index=c("V1"),vSize="V2",type="index",palette=blucol[2:5],title=plottitle)
dev.off()
