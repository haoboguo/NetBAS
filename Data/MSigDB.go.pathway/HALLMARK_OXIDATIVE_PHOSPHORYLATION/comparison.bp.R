####
##This script lists DAVID and NetBAS result side-by-side
##p-values in DAVID are converted to z-scores, and scaled to match the maximum of NetBAS
rm(list=ls())
library(pracma)
library(gplots)

filename <- read.csv("HMK.bp.new.txt", header=T, sep="\t", stringsAsFactors=F)
go.id <- filename$GO.ID[1:10]
go.term <- filename$GO.Term[1:10]
netbas <- filename$zscore[1:10]

netdav <- read.csv("net.dav.txt", header=T, sep=" ", stringsAsFactors=F)
david <- netdav$DAVID[1:10]

#convert DAVID p-values to z-scores
z.david <- c()
for (i in 1:10) {
    z.val <- abs(david[i] - sqrt(2)*erfcinv(2*david[i]))
    z.david <- rbind(z.david, z.val)
}

scale.factor <- max(netbas) / max(z.david[!is.na(z.david)])

z.david <- round(z.david * scale.factor,3)

zscores <- matrix(cbind(netbas, z.david), ncol=2)

colnames(zscores) <- c("NetPAS", "DAVID")
rownames(zscores) <- go.id

colors=c(seq(0,4.9,length=5), seq(5.1,max(netbas),length=5))
my_palette <- colorRampPalette(c("white", "red2"))(n=9)

p.david <- formatC(david, format="E", digits=1)
note <- cbind(netbas, p.david)

png(filename= "Net.DAV.BP10.png", width=3.5, height=6, res=1200, unit="in")
heatmap.2(zscores, col=my_palette, dendrogram='none', breaks=colors,
          colsep = 1:2, rowsep = 1:10, sepcolor="lightgrey", sepwidth=c(0.02,0.02),
          trace='none', Rowv=F, Colv=F, notecol="black", 
          ylab="Biological Process Terms", xlab="",
          margins=c(1,5.5), key.title=NA, key.xlab=NA, key.ylab=NA,
          scale="none", symbreaks=F, symm=F, symkey=F,
          adjCol=c(0.5,0), adjRow=c(0.05,0), srtRow=45, srtCol=0,
          cexRow=0.9, cexCol=0.9, cellnote=note, main="")
dev.off()
