#

library(gplots) #for heatmap.2
cancer <- read.csv("cancer.gene.csv", header=T, stringsAsFactors = F)

cancer.name <- unique(cancer$NAME[!is.na(cancer$SN)])

#calculate association between a pair of two cancers

cancer.dim <- length(cancer.name)

cci.dat <- read.csv("cci.z.csv", header=F, stringsAsFactors = F)
cci.z <- matrix(unlist(cci.dat), nrow=cancer.dim, ncol=cancer.dim)


colnames(cci.z) <- cancer.name
rownames(cci.z) <- cancer.name

my_palette <- colorRampPalette(c("blue2", "white", "red2"))(n = 20)
colors= c(seq(-10,10, length=21))

png("cci.z.heatmap.png", width=12, height=11, res=600, units="in")
heatmap.2(cci.z, col=my_palette, trace='none', breaks=colors, 
          key.xlab=NA, key.title="Z-score", key.ylab=NA, key.xtickfun = NULL, key.ytickfun = NULL,
          srtCol=45, adjCol=c(1,0), dendrogram = "both",
          margins=c(14,18.5), sepwidth=c(0.01,0.01), symbreaks = TRUE,
          sepcolor="grey", colsep=1:length(cancer.name), rowsep=1:length(cancer.name))
dev.off()

