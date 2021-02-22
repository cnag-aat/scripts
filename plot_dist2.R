#!/apps/R/2.14.2/bin/Rscript
args <- commandArgs(trailingOnly = TRUE)
filename <-args[1]
mytitle <-args[2]
lower <-as.numeric(args[3])
upper <-as.numeric(args[4])
pdf(paste(mytitle,".pdf",sep=""))
RF<-read.table(filename,header=FALSE);d<-density(RF[,1],n=5000,from=lower,to=upper);plot(d,xlab="mapping distance (bp)",xlim=c(lower,upper),main=mytitle);abline(v=d$x[which.max(d$y)],col="red");text(round(d$x[which.max(d$y)]), max(d$y), round(d$x[which.max(d$y)]), pos=4);
dev.off()

