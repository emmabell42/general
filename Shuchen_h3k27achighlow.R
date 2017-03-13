#
## Read in all annotated bed files of peaks
#
nice R
library(GenomicRanges)
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
reading.gr <- extend(reading.gr,100,100)
mcols(reading.gr,use.names=T) <- cbind(reading[,1],reading[,6],reading[,8],reading[,16])
assign(toName[i],reading)
assign(paste0(toName[i],".gr"),reading.gr)
}
#
## Determine the median hESC and hNPC H3K27ac Sox2 coverage
#
quantile(as.numeric(as.character(levels(mcols(GSE62193_hESC_H3K27ac_peaks.gr)[,2]))[mcols(GSE62193_hESC_H3K27ac_peaks.gr)[,2]]))
quantile(as.numeric(as.character(levels(mcols(GSE62193_hNPC_H3K27ac_peaks.gr)[,2]))[mcols(GSE62193_hNPC_H3K27ac_peaks.gr)[,2]]))
#GSE62193_hESC_H3K27ac_peaks.gr = 23
#GSE62193_hNPC_H3K27ac_peaks.gr = 30
#
hesc.h3k27ac.high <- GSE62193_hESC_H3K27ac_peaks.gr[which(as.numeric(as.character(mcols(GSE62193_hESC_H3K27ac_peaks.gr)[,2]))>=23)]
hnpc.h3k27ac.high <- GSE62193_hNPC_H3K27ac_peaks.gr[which(as.numeric(as.character(mcols(GSE62193_hNPC_H3K27ac_peaks.gr)[,2]))>=30)]
hesc.h3k27ac.low <- findOverlaps(GSE69479_hESC_Sox2_peaks.gr,hesc.h3k27ac.high)
hesc.h3k37ac.all <- 1:length(GSE62193_hESC_H3K27ac_peaks.gr)
hesc.h3k27ac.low <- hesc.h3k37ac.all[which(!hesc.h3k27ac.all %in% hesc.h3k27ac.low)]
hnpc.h3k27ac.low <- findOverlaps(GSE69479_hnpc_Sox2_peaks.gr,hnpc.h3k27ac.high)
hnpc.h3k37ac.all <- 1:length(GSE62193_hNPC_H3K27ac_peaks.gr)
hnpc.h3k27ac.low <- hnpc.h3k37ac.all[which(!hnpc.h3k27ac.all %in% hnpc.h3k27ac.low)]


hesc.sox2.h3k27ac.high <- subsetByOverlaps(GSE69479_hESC_Sox2_peaks.gr,hnpc.h3k27ac.high)
hesc.sox2.h3k27ac.low <- GSE69479_hESC_Sox2_peaks.gr[!as.numeric(findOverlaps(GSE69479_hESC_Sox2_peaks.gr,hnpc.h3k27ac.low)@queryHits)]
hnpc.sox2.h3k27ac.high <- subsetByOverlaps(GSE69479_hNPC_Sox2_peaks.gr,hesc.h3k27ac.high)
hnpc.sox2.h3k27ac.low <- subsetByOverlaps(GSE69479_hNPC_Sox2_peaks.gr,hesc.h3k27ac.low)
