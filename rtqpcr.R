setwd("\\\\store.ic.ac.uk\\IC\\fom\\surgeryandcancer\\epigenetics-and-development\\Current lab members\\Emma Bell\\Data\\Analysis\\Super Enhancers\\Expression\\")
dct <- read.table("RT-qPCR.txt",sep="\t",stringsAsFactors=F,head=T,row.names=1)
dct <- dct[1:26,1:14]
row.names(dct) <- dct[,1]
dct <- dct[,-1]

ddct <- as.matrix(read.table("ddcts.txt",sep="\t",stringsAsFactors=F,head=T,row.names=1))

cm <- cor(dct)

labs <- c("ESC-Esrrb-/-_empty","ESC-Esrrb-/-_WT","ESC-Esrrb-/-_mu","cEpiSC-Esrrb-/-","cEpiSC-Esrrb+/+","ESC-Esrrb-/-_TK","ESC-Esrrb+/+_TK","ESC-Nanog-/-","ESC-Nanog+/+","ESC-Ncoa3-/-","ESC-Ncoa3+/+","ESC-Esrrb-/-_EB","ESC-Esrrb+/+_EB")
cols <- c("darkgrey","lightgrey","lightgrey","pink","pink","pink","red","lightblue","blue","lightgreen","green","lightyellow","gold")
par(mar=c(10,5,2,2))
for(i in 1:nrow(ddct)){
toName <- paste0(rownames(ddct)[i],"_RelExp.png")
png(toName)
par(mar=c(10,5,5,2))
barplot(ddct[i,],las=2,main=rownames(ddct)[i],names=labs,ylab="Relative expression",col=cols)
dev.off()
}

rescue.ddcts <- dct[,1:3]
rescue.ddcts <- rescue.ddcts/rescue.ddcts[,2]
rescue.ddcts <- na.omit(rescue.ddcts)

# If the ESC-Esrrb-/- mu rescue are more differentiated than the other rescue cells they should be more similar to the cEpiSCs and the empty should be more similar to the ESC-Esrrb-/-.
ddct <- as.matrix(cbind(ddct,rescue.ddcts))
toPlot <- c(4,5,6,8,10,12,14,16)
log.ddct <- log2(ddct)
log.ddct[,toPlot][,1][which(log.ddct[,toPlot][,1]==-Inf)] <- NA
log.ddct[,toPlot][,2][which(log.ddct[,toPlot][,2]==-Inf)] <- NA
cm <- cor(log.ddct[,toPlot],use="pairwise.complete.obs",method="spearman")
dendro <- hclust(as.dist(1-cm))
library(gplots)
heatmap.2(cm,col=bluered(100),trace="none",mar=c(9,9),cexRow=0.8,cexCol=0.8,Rowv=as.dendrogram(dendro),Colv=as.dendrogram(dendro))
