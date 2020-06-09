#
library(gplots) #for heatmap.2

list <- c(1:50)
random.dim <- 50

rri.dat <- read.csv("rri.z.csv", header=F, stringsAsFactors = F)
rri.z <- matrix(unlist(rri.dat), nrow=random.dim, ncol=random.dim)

colnames(rri.z) <- list
rownames(rri.z) <- list

my_palette <- colorRampPalette(c("blue2", "white", "red2"))(n = 20)
#colors= c(seq(-5,-1.5, length=10), seq(-1.49, 1.49, length=10), seq(1.5,5, length=10))
colors = c(seq(-10,10,length=21))

png("rri.z.heatmap.png", width=12, height=11, res=600, units="in")
heatmap.2(rri.z, col=my_palette, trace='none', breaks=colors, 
          key.xlab=NA, key.title="Z-score", key.ylab=NA, key.xtickfun = NULL, key.ytickfun = NULL,
          #srtCol=45, adjCol=c(1,0), dendrogram = "both",
          margins=c(14.5,18), sepwidth=c(0.01,0.01), symbreaks = TRUE,
          sepcolor="grey", colsep=1:random.dim, rowsep=1:random.dim)
dev.off()
