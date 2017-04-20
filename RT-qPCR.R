setwd("\\\\store.ic.ac.uk\\IC\\fom\\surgeryandcancer\\epigenetics-and-development\\Current lab members\\Emma Bell\\Data\\Analysis\\Super Enhancers\\Expression\\")
ddct <- read.table("ddcts.txt",sep="\t",stringsAsFactors=F,head=F,row.names=1)
colnames(ddct) <- ddct[1,]
ddct <- ddct[2:nrow(ddct),]
logfcs <- read.table("logfcs.txt",sep="\t",stringsAsFactors=F,head=F,row.names=1)
colnames(logfcs) <- c("Target",colnames(ddct),"Col")
sems <- read.table("sems.txt",sep="\t",stringsAsFactors=F,head=F,row.names=1)
colnames(sems) <- colnames(ddct)
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
stop("vectors must be same length")
arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

se <- c("lightgreen","pink")

x <- barplot(toplot[grep(paste0(se,collapse="|"),toplot[,17]),])
for(i in c(2,4,7,8,9,11,13,15)){
png(paste0(colnames(logfcs)[i],"_se_log2FC.png"))
toplot <- logfcs[order(logfcs[,i]),]
toplot <- toplot[grep(paste0(se,collapse="|"),toplot[,17]),]
semstoplot <- as.numeric(sems[,i-1])
semstoplot <- semstoplot[order(logfcs[,i])]
semstoplot <- semstoplot[grep(paste0(se,collapse="|"),toplot[,17])]
semstoplot[which(semstoplot==0)] <- NA
semstoplot[is.na(semstoplot)] <- 0.001
ylim <- c(min(na.omit(as.numeric(toplot[,i])-as.numeric(semstoplot))),max(na.omit(as.numeric(toplot[,i])+as.numeric(semstoplot))))*1.1
barplot(toplot[,i],las=2,ylab="Log2(FC)",main=colnames(logfcs)[i],ylim=ylim,names=toplot[,1],col=toplot[,17],cex.names=0.8)
y <- as.numeric(toplot[,i])
y[is.na(y)] <- 0
arrows(as.numeric(x),y-semstoplot,x,y+semstoplot,length=0.05, angle=90, code=3)
grid(nx=NA,ny=NULL)
abline(v=c(6.1,12.1,18.1,24.1),col = "lightgray", lty = "dotted")
dev.off()
}

for(i in c(2,4,7,8,9,11,13,15)){
png(paste0(colnames(logfcs)[i],"_boxplot_log2FC.png"),width=240)
toplot <- logfcs[order(logfcs[,i]),]
toplot <- toplot[grep(paste0(se,collapse="|"),toplot[,17]),]
boxplot(toplot[,i]~toplot[,17],ylab="Log2(FC)",main=colnames(logfcs)[i],names=c("Maintained","Silenced"),col=c("lightgreen","pink"),cex.names=0.8)
dev.off()
}
