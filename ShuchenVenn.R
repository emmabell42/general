#
## Overlapping peaks
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
toSubset <- c("GSE24447_hESC_H3K27me3_peaks.gr","GSE24447_hESC_H3K4me1_peaks.gr","GSE24447_hESC_H3K4me3_peaks.gr")
for(i in 1:3){
subsetting <- subsetByOverlaps(get(toSubset[i]),GSE24447_hESC_H3K27ac_peaks.gr)
assign(gsub(".gr",".ss",toSubset[i]),subsetting)
}
library(ChIPpeakAnno)
library(limma)
res <- makeVennDiagram(Peaks=list(GSE24447_hESC_H3K27me3_peaks.ss,GSE24447_hESC_H3K4me1_peaks.ss,GSE24447_hESC_H3K4me3_peaks.ss),NameOfPeaks=c("H3K27me3", "H3K4me1","H3K4me3"))
png("ESC_histonemods.png")
vennDiagram(res[[2]])
dev.off()
#
## How many enhancers overlap with Sox2 and Oct4/Sox2 and Pax6?
#
h3k27ac <- array(NA,dim=c(length(GSE24447_hESC_H3K27ac_peaks.gr),3))
colnames(h3k27ac) <- toSubset
for(i in 1:ncol(h3k27ac)){
h3k27ac[,i] <- countOverlaps(GSE24447_hESC_H3K27ac_peaks.gr,get(toSubset[i]))
h3k27ac[which(h3k27ac[,i]>1),i] <- 1
}
enhancers <- which(h3k27ac[,1]==0 & h3k27ac[,2]==1 & h3k27ac[,3]==0)
h3k27ac.enhancers <- GSE24447_hESC_H3K27ac_peaks.gr[enhancers] #23507
table(mcols(h3k27ac.enhancers[]))
#
## How many enhancers overlap with Sox2 and Oct4?
#
res <- makeVennDiagram(Peaks=list(h3k27ac.enhancers,GSE69479_hESC_Sox2_peaks.gr,GSE69646_hESC_Oct4_peaks.gr),NameOfPeaks=c("Enhancers","Sox2","Oct4"))
library(limma)
png("hESC_enhancers_Sox2_Oct4_venn.png")
vennDiagram(res[[2]])
dev.off()

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

toExtend <- c("GSE69479_hESC_Sox2_peaks.gr","GSE69646_hESC_Oct4_peaks.gr")
for(i in 1:length(extend)){
gr <- get(toExtend[i])
extended.gr <- extend(gr,100,100)
assign(paste0(toExtend[i],".wide"),extended.gr)
}
res <- makeVennDiagram(Peaks=list(h3k27ac.enhancers,GSE69479_hESC_Sox2_peaks.gr.wide,GSE69646_hESC_Oct4_peaks.gr.wide),NameOfPeaks=c("Enhancers","Sox2","Oct4"))
png("hESC_enhancers_Sox2_Oct4_wide_venn.png")
vennDiagram(res[[2]])
dev.off()
