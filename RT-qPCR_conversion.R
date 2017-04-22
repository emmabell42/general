setwd("\\\\store.ic.ac.uk\\IC\\fom\\surgeryandcancer\\epigenetics-and-development\\Current lab members\\Emma Bell\\Lab work\\qPCR data\\20150826 AJ conversion data\\Final data")
logfcs <- read.table("logfcs.txt",sep="\t",stringsAsFactors=F,head=T)
relexp <- read.table("RelExp.txt",sep="\t",stringsAsFactors=F,head=T)
sems <- read.table("sems.txt",sep="\t",stringsAsFactors=F,head=T)

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
stop("vectors must be same length")
arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}
se <- c("lightgreen","pink")
x <- c(barplot(toplot,beside=T))
for(i in c(5:28)){
png(paste0(colnames(relexp)[i],"_barplot_log2FC.png"))
toplot <- matrix(relexp[c(1,7,2,8,3,9,4,10,5,11,6,12),i],nrow=2,byrow=F)
semstoplot <- as.numeric(sems[,i])
semstoplot[is.na(semstoplot)] <- 0.001
y <- as.numeric(toplot)
y[is.na(y)] <- 0
ylim <- c(0,max(na.omit(as.numeric(toplot)+as.numeric(semstoplot))))*1.1
bp <- barplot(toplot,main=colnames(relexp)[i],ylim=ylim,col=c("lightgrey","darkgrey"),cex.names=0.8,beside=T,xaxt = 'n',yaxt = 'n',xlab="Days",ylab="Relative expression")
grid(nx=NA,ny=NULL)
barplot(toplot,main=colnames(relexp)[i],ylim=ylim,col=c("lightgrey","darkgrey"),cex.names=0.8,beside=T,xaxt = 'n',xlab = '',yaxt = 'n',add=T)
axis(1, at=c(2,5,8,11,14,17), labels=c(0,1,2,3,7,"cEpiSC"),tick=FALSE, las=1, cex.axis=1)
axis(2, at=axTicks(2), cex.axis=0.9, las=2)
arrows(as.numeric(x),y-semstoplot,x,y+semstoplot,length=0.05, angle=90, code=3)
legend("topleft",c("2i","Serum"),fill=c("lightgray","darkgrey"))
dev.off()
}

logfcs.se <- read.table("logfcs_se.txt",sep="\t",head=T)
logfcs.se <- na.omit(logfcs.se)
bp <- boxplot(LogFC~Media*Day,data=logfcs.se,las=2,col=c("lightgrey","darkgrey"),names.arg=NULL,xaxt = 'n',yaxt = 'n',xlab="Days",ylab="Log(FC)")
grid(nx=NA,ny=NULL)
boxplot(LogFC~Media*Day,data=logfcs.se,las=2,col=c("lightgrey","darkgrey"),main="Differential expression of SE associated genes",names.arg=NULL,xaxt = 'n',yaxt = 'n',xlab="Days",ylab="Log(FC)",add=T)
axis(1, at=c(1.5,3.5,5.5,7.5,9.5), labels=c(1,2,3,7,"cEpiSC"),tick=FALSE, las=1, cex.axis=1)
axis(2, at=axTicks(2), cex.axis=0.9, las=2)
legend("topleft",c("2i","Serum"),fill=c("lightgray","darkgrey"),bg="white")
abline(v=c(2.5,4.5,6.5,8.5),lty=3,col="lightgrey")

bp <- boxplot(LogFC~Day,data=logfcs.se[which(logfcs.se$Media=="2i"),],las=2,col="lightgrey",names.arg=NULL,xaxt = 'n',yaxt = 'n',xlab="Days",ylab="Log(FC)")
grid(nx=NA,ny=NULL)
boxplot(LogFC~Day,data=logfcs.se[which(logfcs.se$Media=="2i"),],las=2,col="lightgrey",main="2i",names.arg=NULL,xaxt = 'n',yaxt = 'n',xlab="Days",ylab="Log(FC)",add=T)
axis(1, at=c(1,2,3,4,5), labels=c(1,2,3,7,"cEpiSC"),tick=FALSE, las=1, cex.axis=1)
axis(2, at=axTicks(2), cex.axis=0.9, las=2)

bp <- boxplot(LogFC~Day,data=logfcs.se[which(logfcs.se$Media=="Serum"),],las=2,col="darkgrey",names.arg=NULL,xaxt = 'n',yaxt = 'n',xlab="Days",ylab="Log(FC)")
grid(nx=NA,ny=NULL)
boxplot(LogFC~Day,data=logfcs.se[which(logfcs.se$Media=="Serum"),],las=2,col="darkgrey",main="Serum",names.arg=NULL,xaxt = 'n',yaxt = 'n',xlab="Days",ylab="Log(FC)",add=T)
axis(1, at=c(1,2,3,4,5), labels=c(1,2,3,7,"cEpiSC"),tick=FALSE, las=1, cex.axis=1)
axis(2, at=axTicks(2), cex.axis=0.9, las=2)

bp <- boxplot(LogFC~Setype*Day,data=logfcs.se[which(logfcs.se$Media=="2i"),],las=2,col="lightgrey",names.arg=NULL,xaxt = 'n',yaxt = 'n',xlab="Days",ylab="Log(FC)")
grid(nx=NA,ny=NULL)
boxplot(LogFC~Setype*Day,data=logfcs.se[which(logfcs.se$Media=="2i"),],las=2,col=c("lightgreen","pink"),main="2i",names.arg=NULL,xaxt = 'n',yaxt = 'n',xlab="Days",ylab="Log(FC)",add=T)
axis(1, at=c(1.5,3.5,5.5,7.5,9.5), labels=c(1,2,3,7,"cEpiSC"),tick=FALSE, las=1, cex.axis=1)
axis(2, at=axTicks(2), cex.axis=0.9, las=2)
abline(v=c(2.5,4.5,6.5,8.5),lty=3,col="lightgrey")

bp <- boxplot(LogFC~Setype*Day,data=logfcs.se[which(logfcs.se$Media=="Serum"),],las=2,col="lightgrey",names.arg=NULL,xaxt = 'n',yaxt = 'n',xlab="Days",ylab="Log(FC)")
grid(nx=NA,ny=NULL)
boxplot(LogFC~Setype*Day,data=logfcs.se[which(logfcs.se$Media=="Serum"),],las=2,col=c("lightgreen","pink"),main="Serum",names.arg=NULL,xaxt = 'n',yaxt = 'n',xlab="Days",ylab="Log(FC)",add=T)
axis(1, at=c(1.5,3.5,5.5,7.5,9.5), labels=c(1,2,3,7,"cEpiSC"),tick=FALSE, las=1, cex.axis=1)
axis(2, at=axTicks(2), cex.axis=0.9, las=2)
abline(v=c(2.5,4.5,6.5,8.5),lty=3,col="lightgrey")

#twoi <- matrix(logfcs.se[which(logfcs.se$Media=="2i"),1],ncol=5)
twoi <- t(relexp[which(relexp[,2]=="2i"),5:28])
colnames(twoi) <- c(0,1,2,3,7,"cEpiSC")
#twoi <- 2^twoi
cm <- cor(t(twoi))
dendro <- hclust(as.dist(1-cm))
heatmap.2(twoi,scale="row",Colv=NULL,trace="none",col=bluered(100),Rowv=as.dendrogram(dendro),RowSideColors=as.character(setype),mar=c(7,5))

#serum <- matrix(logfcs.se[which(logfcs.se$Media=="Serum"),1],ncol=5)
serum <- t(relexp[which(relexp[,2]=="Serum"),5:28])
colnames(serum) <- c(0,1,2,3,7,"cEpiSC")
#serum <- 2^serum
cm <- cor(t(serum))
dendro <- hclust(as.dist(1-cm))
heatmap.2(serum,scale="row",Colv=NULL,trace="none",col=bluered(100),Rowv=as.dendrogram(dendro),RowSideColors=as.character(setype),mar=c(7,5))

setype <- c("grey","grey","pink","grey","pink","lightgreen","lightgreen","grey","lightgreen","lightgreen","lightgreen","pink","lightgreen","lightgreen","lightgreen","pink","pink","pink","lightgreen","pink","pink","lightgreen","lightgreen","grey")

silenced <- c("Esrrb","Tbx3","Tfcp2l1","Tet2","Klf2","Klf4","Klf5","Tdh")
maintained <- c("Nanog","Oct4","Klf13","Med13l","Smadcard1","Tet1","Otx2","Lefty1")
ctrl <- c("Dnmt3a","Dnmt3b","Fgf5","Ncoa3","T")

twoi <- twoi[which(rownames(twoi) %in% c(silenced,maintained,ctrl)),]
serum <- serum[which(rownames(serum) %in% c(silenced,maintained,ctrl)),]

logfcs.se <- logfcs.se[which(logfcs.se$Gene %in% c(silenced,maintained,ctrl)),]

setype <- c("grey","grey","pink","grey","pink","lightgreen","lightgreen","grey","lightgreen","pink","pink","lightgreen","pink","pink","lightgreen","lightgreen","grey")
