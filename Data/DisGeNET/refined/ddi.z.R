#
library(gplots) #for heatmap.2


list.file <- read.csv("all.ids.csv", header=T, stringsAsFactors = F)
list <- list.file$ID[which(list.file$GN > 14)]

disorder.dim <- length(list)

ddi.dat <- read.csv("ddi.z.1k.csv", header=F, stringsAsFactors = F)
ddi.z <- matrix(unlist(ddi.dat), nrow=disorder.dim, ncol=disorder.dim)

colnames(ddi.z) <- list
rownames(ddi.z) <- list

my_palette <- colorRampPalette(c("blue2", "white", "red2"))(n = 20)
#colors= c(seq(-5,-1.5, length=10), seq(-1.49, 1.49, length=10), seq(1.5,5, length=10))
colors = c(seq(-10,10,length=21))

png("ddi.z1k.heatmap.png", width=12, height=11, res=1200, units="in")
heatmap.2(ddi.z, col=my_palette, trace='none', breaks=colors, 
          key.xlab=NA, key.title="Z-score", key.ylab=NA, key.xtickfun = NULL, key.ytickfun = NULL,
          srtCol=45, adjCol=c(1,0), dendrogram = "both",
          margins=c(14,18.5)) #, sepwidth=c(0.01,0.01), symbreaks = TRUE,
          #sepcolor="grey", colsep=1:disorder.dim, rowsep=1:disorder.dim)
dev.off()

hm <- heatmap.2(ddi.z, col=my_palette, trace='none', breaks=colors, 
                key.xlab=NA, key.title="Z-score", key.ylab=NA, key.xtickfun = NULL, key.ytickfun = NULL,
                srtCol=45, adjCol=c(1,0), dendrogram = "both",
                margins=c(14,18.5))
hc <- as.hclust(hm$rowDendrogram)
pdf("disorder.tree.pdf", height=5, width=80, paper='special')
plot(hc, xlab="Disorder ID", cex=.8)
dev.off()