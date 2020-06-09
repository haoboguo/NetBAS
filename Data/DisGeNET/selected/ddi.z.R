#
library(gplots) #for heatmap.2

sel <- read.csv("DisGeNET.selected.csv", header=T, stringsAsFactors = F)
list <- unique(sel$ID)

disorder.dim <- 50

ddi.dat <- read.csv("ddi.z.csv", header=F, stringsAsFactors = F)
ddi.z <- matrix(unlist(ddi.dat), nrow=disorder.dim, ncol=disorder.dim)

colnames(ddi.z) <- list
rownames(ddi.z) <- list

my_palette <- colorRampPalette(c("blue2", "white", "red2"))(n = 20)
#colors= c(seq(-5,-1.5, length=10), seq(-1.49, 1.49, length=10), seq(1.5,5, length=10))
colors = c(seq(-10,10,length=21))

png("ddi.z.heatmap.png", width=12, height=11, res=600, units="in")
heatmap.2(ddi.z, col=my_palette, trace='none', breaks=colors, 
          key.xlab=NA, key.title="Z-score", key.ylab=NA, key.xtickfun = NULL, key.ytickfun = NULL,
          srtCol=45, adjCol=c(1,0), dendrogram = "both",
          margins=c(14,18.5), sepwidth=c(0.01,0.01), symbreaks = TRUE,
          sepcolor="grey", colsep=1:disorder.dim, rowsep=1:disorder.dim)
dev.off()
