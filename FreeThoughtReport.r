#Install the R packages you will require.
#Either run R as an administrator to save the packages to your C drive or save them to your personal folder. The latter will create a folder under your My Documents directory to store new R packages).
source("http://bioconductor.org/biocLite.R")
biocLite(c("maps","mapdata","maptools","rgeos","ggplot2","gplots"))n
#Load in the R packages you will need
library(gplots)
library(ggplot2)
library(rgeos)
library(maps)
library(mapdata)
library(maptools)
#Set your working directory to the location your data is stored.
setwd("\\icnas2.cc.ic.ac.uk\\eb2012\\Free Thought Report 2016")
#List the files in your working directory.
list.files()
#There is a single file called "CountryEntryNotes2015.tsv". Read this file into R.
countries <- read.table(file="fot-report-signals-20161125075353.csv",head=T,sep=",",stringsAsFactors=F)
colnames(countries) <- c("Region","Country","Map name","Gov","Edu","Soc","Exp","Summary score")
countries <- countries[order(countries$Map.name),]
http://www.naturalearthdata.com/downloads/
#Read in the .shp files
toRead <- list.files(recursive=T,pattern=".shp")
toName <- c("m10","m110","m50")
for(i in 1:length(toRead)){
tmp <- readShapePoly(toRead[i])
names(tmp) <- toupper(names(tmp))
tmp.f <- fortify(tmp, region = "NAME_LONG")
assign(toName[i],as.data.frame(tmp.f))
}
ids <- unique(tmp.f$id)
incl <- ids[which(ids %in% countries[,2])]
notIncl <- ids[!which(ids %in% countries[,2])]
dataSummary <- as.matrix(cbind(Countries=ids[which(ids %in% countries[,3])],Number=as.numeric(NA)))
countries.f <- merge(tmp.f, countries, by.x = "id", by.y = "Map name")
countries.f <- countries.f[order(countries.f$order),]
qplot(long, lat,data=countries.f, group = group, fill = Gov,geom="polygon")
cbPalette <- c("#adadad","#b6d8ab","#fde69f","#f4b272","#de6468","#440000")
Map <- ggplot(countries.f, aes(long, lat, group = group, fill = Gov)) + geom_polygon(colour="black") + 
scale_color_manual(values=cbPalette) +
coord_equal() + labs(fill = "Status") + 
   ggtitle("Status")
