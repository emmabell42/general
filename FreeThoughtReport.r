#Install the R packages you will require.
#Either run R as an administrator to save the packages to your C drive or save them to your personal folder. The latter will create a folder under your My Documents directory to store new R packages).
source("http://bioconductor.org/biocLite.R")
biocLite(c("maps","mapdata","maptools","rgeos","ggplot2","gplots"))
#Load in the R packages you will need
library(gplots)
library(ggplot2)
library(rgeos)
library(maps)
library(mapdata)
library(maptools)
#Set your working directory to the location your data is stored.
setwd("\\\\icnas2.cc.ic.ac.uk\\eb2012\\Free Thought Report 2016")
#List the files in your working directory.
list.files()
#There is a single file called "CountryEntryNotes2015.tsv". Read this file into R.
countries <- read.table(file="fot-report-signals-20161125075353.csv",head=T,sep=",",stringsAsFactors=F,fill=T)
countries <- countries[order(countries$Map.Name),]
countries$Map.Name[1] <- "Cote dIvoire"
countries <- countries[order(countries$Map.Name),]
countries$Const.Govt[which(is.na(countries$Const.Govt))] <- 0
countries$Edu.Child[which(is.na(countries$Edu.Child))] <- 0
countries$Society.Comm[which(is.na(countries$Society.Comm))] <- 0
countries$Expression[which(is.na(countries$Expression))] <- 0
#http://www.naturalearthdata.com/downloads/
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
tmp.f$id[which(tmp.f$id=="Côte d'Ivoire")] <- "Cote dIvoire"
incl <- ids[which(ids %in% countries[,2])]
notIncl <- ids[!which(ids %in% countries[,2])]
toMap <- paste("map",toName,sep="_")
for(i in 1:3){
   map <- get(toName[i])
   tmp <- merge(map, countries, by.x = "id", by.y = "Map.Name")
   tmp <- tmp[order(tmp$order),]
   assign(toMap[i],tmp)
}
#qplot(long, lat,data=countries.f, group = group, fill = Const.Govt,geom="polygon")
cbPalette <- c("#adadad","#b6d8ab","#fde69f","#f4b272","#de6468","#440000")
   png(paste0("map_m110","_Const.Govt.png"))
   ggplot(map_m110, aes(long, lat, group = group, fill = factor(Const.Govt))) + geom_polygon(colour="darkgrey") + 
   scale_fill_manual(values=cbPalette) +
   coord_equal() + labs(fill = "Status") + 
   ggtitle("Const.Govt")
   dev.off()
   png(paste0(toMap[i],"_Edu.Child.png"))
   ggplot(map, aes(long, lat, group = group, fill = factor(Edu.Child))) + geom_polygon(colour="darkgrey") + 
   scale_fill_manual(values=cbPalette) +
   coord_equal() + labs(fill = "Status") + 
   ggtitle("Edu.Child")
   png(paste0(toMap[i],"_Society.Comm.png"))
   ggplot(map, aes(long, lat, group = group, fill = factor(Society.Comm))) + geom_polygon(colour="darkgrey") + 
   scale_fill_manual(values=cbPalette) +
   coord_equal() + labs(fill = "Status") + 
   ggtitle("Society.Comm")
   png(paste0(toMap[i],"_Expression.png"))
   ggplot(map, aes(long, lat, group = group, fill = factor(Expression))) + geom_polygon(colour="darkgrey") + 
   scale_fill_manual(values=cbPalette) +
   coord_equal() + labs(fill = "Status") + 
   ggtitle("Expression")
   png(paste0(toMap[i],"_Summary.png"))
   ggplot(map, aes(long, lat, group = group, fill = factor(Summary.Score))) + geom_polygon(colour="darkgrey") + 
   coord_equal() + labs(fill = "Status") + 
   ggtitle("Summary.Score")
   dev.off()
}
