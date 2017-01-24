nice R
library(GenomicRanges)
library(ShortRead)
#
## Create list of files to read in
#
bedName <- list.files(recursive=T)[grep(".narrowPeak",list.files(recursive=T))] # Will find narrowPeak files in the current and sub-directories
path <- paste0(getwd(),toRead)
grName <- paste0(toRead,".gr")
#
## Read in bed files
#
for(i in 1:length(bedName)){
reading <- read.table(path[i],head=F,stringsAsFactors=F,comment.char="",quote="")
reading.gr <- with(reading,GRanges(reading[,1],IRanges(reading[,2],reading[,3])))
assign(bedName[i],reading)
assign(grName[i],reading.gr)
}
#
## Create table of promoters
#
library(biomaRt)
ensembl <- useMart(host="feb2014.archive.ensembl.org",biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")
geneInfo <- getBM(attributes=c("hgnc_symbol","chromosome_name","start_position","end_position","strand"),mart=ensembl)
geneInfo <- geneInfo[which(!geneInfo$hgnc_symbol==""),]
chromosomes <- c(1:22,"X","Y")
geneInfo <- geneInfo[which(geneInfo$chromosome_name %in% chromosomes),]
geneInfo <- geneInfo[!duplicated(geneInfo$hgnc_symbol),]
geneInfo.TSS <- geneInfo # Define your promoter regions however you like on lines 31 and 33
geneInfo.TSS[which(geneInfo.TSS$strand==1),"end_position"] <- geneInfo.TSS[which(geneInfo.TSS$strand==1),"start_position"]
geneInfo.TSS[which(geneInfo.TSS$strand==1),"start_position"] <- geneInfo.TSS[which(geneInfo.TSS$strand==1),"start_position"]-2000
geneInfo.TSS[which(geneInfo.TSS$strand==-1),"start_position"] <- geneInfo.TSS[which(geneInfo.TSS$strand==-1),"end_position"]
geneInfo.TSS[which(geneInfo.TSS$strand==-1),"end_position"] <- geneInfo.TSS[which(geneInfo.TSS$strand==-1),"end_position"]+2000
promoters.gr <- with(geneInfo.TSS,GRanges(chromosome_name,IRanges(start_position,end_position)))
proms <- array(NA,dim=c(nrow(geneInfo.TSS),4))
colnames(proms) <- # Histone modifications and conditions
#
## Overlap promoters with bed files
#
gr <- c() # List GRanges objects in the order you want the columns in the proms to be filled
for(i in 1:3){
proms[,i] <- countOverlaps(promoters.gr.wide,get(gr[i]))
}
