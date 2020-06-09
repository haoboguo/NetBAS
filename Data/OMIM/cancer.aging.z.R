#
library(gplots) #for heatmap.2
cancer <- read.csv("cancer.gene.csv", header=T, stringsAsFactors = F)
cancer.name <- unique(cancer$NAME[!is.na(cancer$SN)])

cancer.dim <- length(cancer.name)

cai.dat <- read.csv("cai.z.csv", header=F, stringsAsFactors = F)
cai.z <- matrix(unlist(cai.dat), ncol=cancer.dim, nrow=2)

colnames(cai.z) <- cancer.name
rownames(cai.z) <- c("Over Expression in Aging", "Under Expression in Aging")

my_palette <- colorRampPalette(c("blue2", "white", "red2"))(n = 29)
colors= c(seq(-3,-0.51, length=10), seq(-0.49, 0.49, length=10), seq(0.51,3, length=10))

png("cai.z.heatmap.png", width=10, height=4, res=600, units="in")
heatmap.2(cai.z, col=my_palette, trace='none', breaks=colors, 
          key.xlab=NA, key.title="Z-score", key.ylab=NA,
          keysize=1, symbreaks=TRUE,
          srtCol=45, adjCol=c(1,0), dendrogram = "both",
          margins=c(8,12), sepwidth=c(0.01,0.01), cexRow = 0.8, cexCol = 0.8,
          sepcolor="grey", colsep=1:cancer.dim, rowsep=1:cancer.dim)
dev.off()


