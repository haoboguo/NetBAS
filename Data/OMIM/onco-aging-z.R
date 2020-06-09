library(gplots) #for heatmap.2

onco <- read.csv("onco.signature.uponly.csv", sep=",", header=T, stringsAsFactors = F)
oncolist <- colnames(onco)
cancer.dim <- length(oncolist)

oai.dat <- read.csv("oai.z.csv", header=F, stringsAsFactors = F)
oai.z <- matrix(unlist(oai.dat), ncol=cancer.dim, nrow=2)

colnames(oai.z) <- oncolist
rownames(oai.z) <- c("Over Expression in Aging", "Under Expression in Aging")

my_palette <- colorRampPalette(c("blue2", "white", "red2"))(n = 29)
colors= c(seq(-3,3, length=30))

png("oai.z.heatmap.png", width=12, height=4, res=600, units="in")
heatmap.2(oai.z, col=my_palette, trace='none', breaks=colors, 
          key.xlab=NA, key.title="Z-score", key.ylab=NA,
          keysize=1, symbreaks=TRUE,
          srtCol=45, adjCol=c(1,0), dendrogram = "both",
          margins=c(10,9), sepwidth=c(0.01,0.01), cexRow = 0.8, cexCol = 0.8,
          sepcolor="grey", colsep=1:cancer.dim, rowsep=1:cancer.dim)
dev.off()
