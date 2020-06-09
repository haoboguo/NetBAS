#
library(gplots) #for heatmap.2

onco <- read.csv("onco.signature.uponly.csv", sep=",", header=T, stringsAsFactors = F)
oncolist <- colnames(onco)
cancer.dim <- length(oncolist)

ooi.dat <- read.csv("ooi.z.csv", header=F, stringsAsFactors = F)
ooi.z <- matrix(unlist(ooi.dat), nrow=cancer.dim, ncol=cancer.dim)

colnames(ooi.z) <- oncolist
rownames(ooi.z) <- oncolist

my_palette <- colorRampPalette(c("blue2", "white", "red2"))(n = 20)
#colors= c(seq(-5,-1.5, length=10), seq(-1.49, 1.49, length=10), seq(1.5,5, length=10))
colors = c(seq(-10,10,length=21))

png("ooi.z.heatmap.png", width=12, height=11, res=600, units="in")
heatmap.2(ooi.z, col=my_palette, trace='none', breaks=colors, 
          key.xlab=NA, key.title="Z-score", key.ylab=NA, key.xtickfun = NULL, key.ytickfun = NULL,
          srtCol=45, adjCol=c(1,0), dendrogram = "both",
          margins=c(14.5,18), sepwidth=c(0.01,0.01), symbreaks = TRUE,
          sepcolor="grey", colsep=1:cancer.dim, rowsep=1:cancer.dim)
dev.off()