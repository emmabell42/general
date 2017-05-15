targets <- unique(cnrq$Target)
clones <- unique(cnrq$Clone)

for(i in 1:length(targets)){
target <- targets[i]
toplot <- cnrq[which(cnrq$Target==target),]
average <- sumtab[which(sumtab$Target==target),]
ylim.2i <- c(0,max(toplot$RNCQ[which(toplot$Condition=="2i")],na.rm=T)*1.1)
ylim.serum <- c(0,max(toplot$RNCQ[which(toplot$Condition=="Serum")],na.rm=T)*1.1)
xlim <-c(0,10)
png(paste0(target,"_2i_RelExp.png"))
lp <- plot(RNCQ~Day,data=toplot[which(toplot$Condition=="2i" & toplot$Clone==clones[1]),],xlim=xlim,ylim=ylim.2i,type="o",pch=20,col="grey",main=target,names.arg=NULL,xaxt = 'n',yaxt = 'n',xlab="Day",ylab="Relative expression")
grid(nx=NA,ny=NULL)
points(RNCQ~Day,data=toplot[which(toplot$Condition=="2i" & toplot$Clone==clones[2]),],type="o",pch=20,col="grey",names.arg=NULL,xaxt = 'n',yaxt = 'n')
points(RNCQ~Day,data=toplot[which(toplot$Condition=="2i" & toplot$Clone==clones[3]),],type="o",pch=20,col="grey",names.arg=NULL,xaxt = 'n',yaxt = 'n')
points(RNCQ~Day,data=toplot[which(toplot$Condition=="2i" & toplot$Clone==clones[4]),],type="o",pch=20,col="grey",names.arg=NULL,xaxt = 'n',yaxt = 'n')
points(RNCQ~Day,data=toplot[which(toplot$Condition=="2i" & toplot$Clone==clones[5]),],type="o",pch=20,col="grey",names.arg=NULL,xaxt = 'n',yaxt = 'n')
points(RNCQ~Day,data=toplot[which(toplot$Condition=="2i" & toplot$Clone==clones[6]),],type="o",pch=20,col="grey",names.arg=NULL,xaxt = 'n',yaxt = 'n')
points(RelExp~Day,data=average[which(average$Condition=="2i"),],type="o",pch=16,names.arg=NULL,xaxt = 'n',yaxt = 'n')
axis(1, at=c(0,1,2,3,7,10), labels=c(0,1,2,3,7,"cEpiSC"),tick=T, las=1, cex.axis=1)
axis(2, at=axTicks(2), cex.axis=0.9, las=2)
dev.off()
png(paste0(target,"_Serum_RelExp.png"))
lp <- plot(RNCQ~Day,data=toplot[which(toplot$Condition=="Serum" & toplot$Clone==clones[1]),],xlim=xlim,ylim=ylim.serum,type="o",pch=20,col="grey",main=target,names.arg=NULL,xaxt = 'n',yaxt = 'n',xlab="Day",ylab="Relative expression")
grid(nx=NA,ny=NULL)
points(RNCQ~Day,data=toplot[which(toplot$Condition=="Serum" & toplot$Clone==clones[2]),],type="o",pch=20,col="grey",names.arg=NULL,xaxt = 'n',yaxt = 'n')
points(RNCQ~Day,data=toplot[which(toplot$Condition=="Serum" & toplot$Clone==clones[3]),],type="o",pch=20,col="grey",names.arg=NULL,xaxt = 'n',yaxt = 'n')
points(RNCQ~Day,data=toplot[which(toplot$Condition=="Serum" & toplot$Clone==clones[4]),],type="o",pch=20,col="grey",names.arg=NULL,xaxt = 'n',yaxt = 'n')
points(RNCQ~Day,data=toplot[which(toplot$Condition=="Serum" & toplot$Clone==clones[5]),],type="o",pch=20,col="grey",names.arg=NULL,xaxt = 'n',yaxt = 'n')
points(RNCQ~Day,data=toplot[which(toplot$Condition=="Serum" & toplot$Clone==clones[6]),],type="o",pch=20,col="grey",names.arg=NULL,xaxt = 'n',yaxt = 'n')
points(RelExp~Day,data=average[which(average$Condition=="Serum"),],type="o",pch=16,names.arg=NULL,xaxt = 'n',yaxt = 'n')
axis(1, at=c(0,1,2,3,7,10), labels=c(0,1,2,3,7,"cEpiSC"),tick=T, las=1, cex.axis=1)
axis(2, at=axTicks(2), cex.axis=0.9, las=2)
dev.off()
}

