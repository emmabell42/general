#
## Overlapping peaks
#
#
## Read in all annotated bed files of peaks
#
nice R
library(GenomicRanges)
path <- paste0(getwd(),"/",list.files(recursive=T)[grep("annotation.txt",list.files(recursive=T))])
toName <- gsub("_annotation.txt","",list.files(recursive=T)[grep("annotation.txt",list.files(recursive=T))])
toName <- gsub("MACS2/","",toName)
for(i in 1:length(path)){
reading <- read.table(path[i],sep="\t",head=T,comment.char="",quote="",stringsAsFactors=F)
reading$Annotation[grep("intron",reading$Annotation)] <- "intron"
reading$Annotation[grep("promoter-TSS",reading$Annotation)] <- "promoter-TSS"
reading$Annotation[grep("exon",reading$Annotation)] <- "exon"
reading$Annotation[grep("non-coding",reading$Annotation)] <- "non-coding"
reading$Annotation[grep("TTS",reading$Annotation)] <- "TTS"
cat(toName[i]," ",table(reading$Annotation),"\n")
reading.gr <- with(reading,GRanges(reading[,2],IRanges(reading[,3],reading[,4])))
mcols(reading.gr,use.names=T) <- cbind(reading[,1],reading[,8],reading[,16])
assign(toName[i],reading)
assign(paste0(toName[i],".gr"),reading.gr)
}
#
## For both hESCs and hNPCs, subset H3K27me3, H3K4me1 and H3k4me3 to those peaks that overlap with H3K27ac
#
toSubset <- list(c("GSE24447_hESC_H3K27me3_peaks.gr","GSE24447_hESC_H3K4me1_peaks.gr","GSE24447_hESC_H3K4me3_peaks.gr"),c("GSE24447_hNPC_H3K27me3_peaks.gr","GSE24447_hNPC_H3K4me1_peaks.gr","GSE24447_hNPC_H3K4me3_peaks.gr"))
names(toSubset) <- c("GSE24447_hESC_H3K27ac_peaks.gr","GSE24447_hNPC_H3K27ac_peaks.gr")
for(i in 1:length(toSubset)){
set <- toSubset[[i]]
h3k27ac <- get(names(toSubset)[[i]])
for(j in 1:length(set)){
subsetting <- subsetByOverlaps(get(set[j]),h3k27ac)
assign(gsub(".gr",".ss",set[j]),subsetting)
}
}
#
## Create counts tables of overlapping histone modification peaks
#
library(ChIPpeakAnno)
library(limma)
subsetted <- toSubset
names(subsetted) <- c("hESC_H3K27ac","hNPC_H3K27ac")
res.hist <- as.list(rep(NA,2))
for(i in 1:length(res.hist)){
subsetted[[i]] <- gsub(".gr",".ss",subsetted[[i]])
res.hist[[i]] <- makeVennDiagram(Peaks=list(get(subsetted[[i]][1]),get(subsetted[[i]][2]),get(subsetted[[i]][3])),NameOfPeaks=c("H3K27me3", "H3K4me1","H3K4me3"))
#png(paste0(names(subsetted)[[i]],"_venn.png"))
#vennDiagram(res.hist[[i]][[2]])
#dev.off()
}
#
## Define enhancers in both cell lines
#
h3k27ac <- list(array(NA,dim=c(length(GSE24447_hESC_H3K27ac_peaks.gr),3)),array(NA,dim=c(length(GSE24447_hNPC_H3K27ac_peaks.gr),3)))
enhancers <- as.list(rep(NA,2))
names(enhancers) <- gsub("H3K27ac","enhancers",names(subsetted))
for(i in 1:length(h3k27ac)){
colnames(h3k27ac[[i]]) <- subsetted[[i]]
	for(j in 1:ncol(h3k27ac[[i]])){
	h3k27ac[[i]][,j] <- countOverlaps(get(names(toSubset)[[i]]),get(subsetted[[i]][j]))
	h3k27ac[[i]][which(h3k27ac[[i]][,j]>1),j] <- 1
	}
enh <- which(h3k27ac[[i]][,1]==0 & h3k27ac[[i]][,2]==1 & h3k27ac[[i]][,3]==0)
enhancers[[i]] <- get(names(toSubset)[[i]])[enh]
}
table(mcols(enhancers[[1]])[2])
#> table(mcols(enhancers[[1]])[2])
#
#        exon   Intergenic       intron promoter-TSS          TTS
#        1253         7840        13548          373          489
table(mcols(enhancers[[2]])[2])
#> table(mcols(enhancers[[2]])[2])
#
#        exon   Intergenic       intron promoter-TSS          TTS
#         181         2455         3691           83           96
#
#
#
## How many enhancers overlap with Sox2 and Oct4/Sox2 and Pax6?
#
tfs <- list(c("GSE69479_hESC_Sox2_peaks.gr","GSE69646_hESC_Oct4_peaks.gr"),c("GSE69479_hNPC_Sox2_peaks.gr","Pax6.homer.gr"))
names(tfs) <- c("hESC","hNPC")
tfs.ext <- tfs
extend <- function(x, upstream=0, downstream=0)     
{
    if (any(strand(x) == "*"))
        warning("'*' ranges were treated as '+'")
    on_plus <- strand(x) == "+" | strand(x) == "*"
    new_start <- start(x) - ifelse(on_plus, upstream, downstream)
    new_end <- end(x) + ifelse(on_plus, downstream, upstream)
    ranges(x) <- IRanges(new_start, new_end)
    trim(x)
}
for(i in 1:length(tfs)){
tfs.ext[[i]] <- gsub(".gr",".gr.wide",tfs[[i]])
	for(j in 1:length(tfs[[i]])){
	gr <- get(tfs[[i]][j])
	extended.gr <- extend(gr,100,100)
	assign(paste0(tfs[[i]][j],".wide"),extended.gr)
	}
}
res.enh <- as.list(rep(NA,2))
names(res.enh) <- c("hESC_enhancers_Sox2_Oct4","hNPC_enhancers_Sox2_Pax6")
for(i in 1:length(res.enh)){
res.enh[[i]] <- makeVennDiagram(Peaks=list(enhancers[[i]],get(tfs.ext[[i]][1]),get(tfs.ext[[i]][2])),NameOfPeaks=c("Enhancers",tfs[[i]]))
png(paste0(names(enhancers)[[i]],"_tfs_venn.png"))
vennDiagram(res.enh[[i]][[2]])
dev.off()
}
#
## Define enhancers that overlap TFs in both cell lines
#
enhancers.tf <- list(array(NA,dim=c(length(enhancers[[1]]),2)),array(NA,dim=c(length(enhancers[[2]]),2)))
names(enhancers.tf) <- names(res.enh)
geneLists <- as.list(rep(NA,2))
names(geneLists) <- names(res.enh)
for(i in 1:length(enhancers.tf)){
colnames(enhancers.tf[[i]]) <- tfs.ext[[i]]
	for(j in 1:ncol(enhancers.tf[[i]])){
	enhancers.tf[[i]][,j] <- countOverlaps(enhancers[[i]],get(tfs.ext[[i]][j]))
	enhancers.tf[[i]][which(enhancers.tf[[i]][,j]>1),j] <- 1
	}
enh <- which(enhancers.tf[[i]][,1]==1 & enhancers.tf[[i]][,2]==1)
geneLists[[i]] <- enhancers[[i]][enh]
df <- data.frame(chr=seqnames(geneLists[[i]]),start=start(geneLists[[i]])-1,end=end(geneLists[[i]]),names=mcols(geneLists[[i]])[1],locations=mcols(geneLists[[i]])[2],gene=mcols(geneLists[[i]])[3])
assign(names(geneLists)[[i]],df)
write.table(df,paste0(names(geneLists)[[i]],".txt"),sep="\t",row.names=F,col.names=F,quote=F)
}
#> table(mcols(geneLists[[1]])[2])
#
#        exon   Intergenic       intron promoter-TSS          TTS
#          33          875          805           13           21
#> table(mcols(geneLists[[2]])[2])
#
#        exon   Intergenic       intron promoter-TSS          TTS
#          20          449          524            7            5
#
## All hESC Sox2/Oct4 and Sox2/Pax6 overlaps
#
res.tf <- as.list(rep(NA,2))
names(res.tf) <- c("hESC_Sox2_Oct4","hNPC_Sox2_Pax6")
for(i in 1:length(res.tf)){
res.tf[[i]] <- makeVennDiagram(Peaks=list(get(tfs.ext[[i]][1]),get(tfs.ext[[i]][2])),NameOfPeaks=c(tfs.ext[[i]]))
#png(paste0(names(res.tf)[[i]],"_venn.png"))
#vennDiagram(res.tf[[i]][[2]])
#dev.off()
}
tf.overlaps <- list(rep(NA,length(get(tfs.ext[[1]][1]))),rep(NA,length(get(tfs.ext[[2]][1]))))
names(tf.overlaps) <- names(res.tf)
geneLists <- as.list(rep(NA,2))
names(geneLists) <- names(res.tf)
for(i in 1:length(tf.overlaps)){
tf.overlaps[[i]] <- countOverlaps(get(tfs.ext[[i]][1]),get(tfs.ext[[i]][2]))
tf.overlaps[[i]][which(tf.overlaps[[i]]>1)] <- 1
geneLists[[i]] <- get(tfs.ext[[i]][1])[which(tf.overlaps[[i]]==1)]
df <- data.frame(chr=seqnames(geneLists[[i]]),start=start(geneLists[[i]])-1,end=end(geneLists[[i]]),names=mcols(geneLists[[i]])[1],locations=mcols(geneLists[[i]])[2],gene=mcols(geneLists[[i]])[3])
assign(names(geneLists)[[i]],df)
#write.table(df,paste0(names(geneLists)[[i]],".txt"),sep="\t",row.names=F,col.names=F,quote=F)
}
##############################
#
## How many hESC Sox2-Oct4 co-bound regions overlap with hNPC Sox2 and hNPC Pax6?
#
##############################
library(VennDiagram)
res1 <- makeVennDiagram(Peaks=list(geneLists[[i]],GSE69479_hNPC_Sox2_peaks.gr.wide,Pax6.homer.gr.wide),NameOfPeaks=c("hESC Sox2-Oct4","hNPC Sox2","hNPC Pax6"))
png("venn_hESC_Sox2_Oct4_hNPC_Sox2_Pax6.png")
vennDiagram(res[[2]])
dev.off()
res <- makeVennDiagram(Peaks=list(GSE69479_hESC_Sox2_peaks.gr.wide,GSE69646_hESC_Oct4_peaks.gr.wide,GSE69479_hNPC_Sox2_peaks.gr.wide,Pax6.homer.gr.wide),NameOfPeaks=c("hESC Sox2","hESC Oct4","hNPC Sox2","hNPC Pax6"))

####################################
#
## Heatmaps
#
####################################
setwd("/data/emmabell42/seq/Shuchen/coverage")
hesc.cov <- read.table("20170126_hESC_enhancers",sep="\t",head=T,row.names=1)
hnpc.cov <- read.table("20170126_hNPC_enhancers",sep="\t",head=T,row.names=1)
tags <- list.files("../tagdir")[1:19]
ncol(cov)/161 #19
covList <- as.list(rep(NA,19))


aveCov <- list(array(NA,dim=c(nrow(hesc.cov),19)),array(NA,dim=c(nrow(hnpc.cov),19)))
toCalc <- c("hesc.cov","hnpc.cov")
colnames(aveCov[[1]]) <- tags
colnames(aveCov[[2]]) <- tags
for(i in 1:length(aveCov)){
cov <- get(toCalc[i])
	for(j in 1:nrow(cov)){
	cat(i,"Calculating means for row",j,"\n",sep=" ")
	  for(k in 1:ncol(aveCov[[i]])){
		#cat("Calculating mean of",tags[k],"\n",sep=" ")
		lastCol <- k*41
		firstCol <- lastCol-40
		aveCov[[i]][j,k] <- mean(as.numeric(cov[j,firstCol:lastCol]))
		}
	}
}

library(gplots)
cm <- as.list(rep(NA,2))
names(cm) <- c("hESC_enhancers","hNPC_enhancers")
dendro <- as.list(rep(NA,2))
toSelect <- list("hESC",c("hNPC","Pax6"))
for(i in 1:(length(cm))){
cm[[i]] <- cor(aveCov[[i]][,grep(paste0(toSelect[[i]],collapse="|"),colnames(aveCov[[i]]))],use="pairwise.complete.obs",method="spearman")
dendro[[i]] <- hclust(as.dist(1-cm[[i]]))
png(paste0(names(cm)[i],"_heatmap.png"))
heatmap.2(cm[[i]],Rowv=as.dendrogram(dendro[[i]]),Colv=as.dendrogram(dendro[[i]]),trace="none",col=bluered(100),mar=c(15,15))
dev.off()
}

dendro <- hclust(as.dist(1-cm))
png("hESC_enhancers.png")
heatmap.2(,Rowv=as.dendrogram(dendro),Colv=as.dendrogram(dendro),trace="none",col=bluered(100),mar=c(15,15))
dev.off()

#
## Motif analysis at hESC Sox2-Oct4 and hNPC Sox2-Pax6
#
PATH=$PATH:/data/seqtools/homer/bin/
PATH=$PATH:/data/seqtools/weblogo/
/data/seqtools/homer/bin/findMotifsGenome.pl hESC_Sox2_Oct4.txt hg19 MotifEnrichment_hESC_Sox2_Oct4 -preparsedDir ../../resources/HOMER
/data/seqtools/homer/bin/findMotifsGenome.pl hNPC_Sox2_Pax6.txt hg19 MotifEnrichment_hNPC_Sox2_Pax6 -preparsedDir ../../resources/HOMER
/data/seqtools/homer/bin/findMotifsGenome.pl hESC_Sox2_Oct4.txt hg19 MotifEnrichment_hESC_Sox2_Oct4_vs_hNPC -bg hNPC_Sox2_Pax6.txt
/data/seqtools/homer/bin/findMotifsGenome.pl hNPC_Sox2_Pax6.txt hg19 MotifEnrichment_hNPC_Sox2_Pax6_vs_hESC -bg hESC_Sox2_Oct4.txt

#
## 20170626
#
library(GenomicRanges)
library(ChIPpeakAnno)

# By ChIP-seq peak

toRead <- list.files()
toName <- gsub("_peaks.narrowPeak","",toRead)
toName <- gsub("_peaks.broadPeak","",toName)
for(i in 1:length(toRead)){
tmp <- read.table(toRead[i],head=F,stringsAsFactors=F,comment.char="",quote="",sep="\t")
colnames(tmp)[1:3] <- c("chr","start","end")
tmp.gr <- with(tmp, GRanges(chr, IRanges(start=start, end=end)))
assign(toName[i],tmp)
assign(paste0(toName[i],".gr"),tmp.gr)
}
Pax6.homer.bed$chr <- gsub(" ","",Pax6.homer.bed$chr)
Pax6.homer.bed.gr <- with(Pax6.homer.bed, GRanges(chr, IRanges(start=start, end=end)))
gr <- paste0(toName,".gr")
so <- subsetByOverlaps(GSE69479_hESC_Sox2.gr,GSE69646_hESC_Oct4.gr)
sp <- subsetByOverlaps(GSE69479_hNPC_Sox2.gr,Pax6.homer.bed.gr)
gr <- c("so","sp","GSE62193_hESC_H3K27ac.gr","GSE62193_hNPC_H3K27ac.gr")
venn <- makeVennDiagram(peaks=list(so,sp,GSE62193_hESC_H3K27ac.gr,GSE62193_hNPC_H3K27ac.gr),NameOfPeaks=gr)

# By gene

toRead <- list.files()[grep("annotation",list.files())]
toName <- gsub("_peaks_annotation.txt",".annot",toRead)
genes <- c()
for(i in 1:length(toRead)){
tmp <- read.table(toRead[i],head=T,stringsAsFactors=F,comment.char="",quote="",sep="\t")
genes <- c(genes,tmp$Gene.Name)
tmp.gr <- with(tmp, GRanges(Chr, IRanges(start=Start, end=End),mcols=Gene.Name))
assign(toName[i],tmp)
assign(paste0(toName[i],".gr"),tmp.gr)
}
genes <- unique(genes)
genes <- sort(genes)
genes <- genes[2:length(genes)]
venn.genes <- array(NA,dim=c(length(genes),4))
rownames(venn.genes) <- genes 
colnames(venn.genes) <- gr
so <- subsetByOverlaps(GSE69479_hESC_Sox2.annot.gr,GSE69646_hESC_Oct4.annot.gr)
sp <- subsetByOverlaps(GSE69479_hNPC_Sox2.annot.gr,Pax6.homer_annotation.txt.gr)
gr <- c("so","sp","GSE62193_hESC_H3K27ac.annot.gr","GSE62193_hNPC_H3K27ac.annot.gr")
for(i in 1:4){
regs <- get(gr[i])
venn.genes[,i] <- as.numeric(rownames(venn.genes) %in% mcols(regs)$mcols)
}
tmp <- rowSums(venn.genes)
venn.genes <- venn.genes[which(tmp>0),]
