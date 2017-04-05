setwd("\\\\store.ic.ac.uk\\IC\\fom\\surgeryandcancer\\epigenetics-and-development\\Current lab members\\Emma Bell\\Data\\Analysis\\Super Enhancers\\Expression\\")
ddct <- read.table("ddcts.txt",sep="\t",stringsAsFactors=F,head=F,row.names=1)
colnames(ddct) <- ddct[1,]
ddct <- ddct[2:nrow(ddct),]
logfcs <- read.table("logfc.txt",sep="\t",stringsAsFactors=F,head=F,row.names=1)
colnames(logfcs) <- labs
sems <- read.table("sems.txt",sep="\t",stringsAsFactors=F,head=F,row.names=1)
geomean = function(x, na.rm=TRUE, zero.propagate = FALSE){
  if(any(x < 0, na.rm = TRUE)){
    return(NaN)
  }
  if(zero.propagate){
    if(any(x == 0, na.rm = TRUE)){
      return(0)
    }
    exp(mean(log(x), na.rm = na.rm))
  } else {
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
}
toReplace <- c("_Sample.*","_p.*")
labs <- unique(gsub(paste0(toReplace,collapse="|"),"",colnames(ddct)))
geomeans <- array(NA,dim=c(nrow(ddct),length(labs)))
colnames(geomeans) <- labs
rownames(geomeans) <- rownames(ddct)
sems <- geomeans

for(i in 1:ncol(geomeans)){
toCalc <- labs[i]
	for(j in 1:nrow(geomeans)){
	values <- as.numeric(ddct[j,grep(toCalc,colnames(ddct))])
	geomeans[j,i] <- geomean(values)
	sems[j,i] <- sd(log2(values))/sqrt(length(values))
	}
}

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
stop("vectors must be same length")
arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

arrows(x0, y0, x1 = x0, y1 = y0
x <- barplot(as.numeric(logfcs[1,]))
for(i in 1:nrow(logfcs)){
png(paste0(rownames(logfcs)[i],"_log2FC.png"))
ylim <- c(min(na.omit(as.numeric(logfcs[i,])-as.numeric(sems[i,]))),max(na.omit(as.numeric(logfcs[i,])+as.numeric(sems[i,]))))*1.1
barplot(as.numeric(logfcs[i,]),las=2,ylab="Log2(FC)",main=rownames(logfcs)[i],ylim=ylim,names=labs,cex.names=0.7)
y <- as.numeric(logfcs[i,])
y[is.na(y)] <- 0
semstoplot <- as.numeric(sems[i,])
semstoplot[which(semstoplot==0)] <- NA
semstoplot[is.na(semstoplot)] <- 0.001
arrows(as.numeric(x),y-semstoplot,x,y+semstoplot,length=0.05, angle=90, code=3)
grid(nx=NA,ny=NULL)
abline(v=c(4.9,8.5,10.9,13.3,15.7),col = "lightgray", lty = "dotted")
dev.off()
}
