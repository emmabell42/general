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
png(paste0(names(subsetted)[[i]],"_venn.png"))
vennDiagram(res.hist[[i]][[2]])
dev.off()
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
}


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
dendro <- as.list(rep(NA,2))
for(i in 1:(length(cm))){
cm[[i]] <- cor(aveCov[[i]],use="pairwise.complete.obs")
dendro[[i]] <- hclust(as.dist(1-cm[[i]]))
}

dendro <- hclust(as.dist(1-cm))
png("hESC_enhancers.png")
heatmap.2(,Rowv=as.dendrogram(dendro),Colv=as.dendrogram(dendro),trace="none",col=bluered(100),mar=c(15,15))
dev.off()

#
## Which enhancers overlap with Sox2-Oct4 and Sox2-Pax6? 
#

	
