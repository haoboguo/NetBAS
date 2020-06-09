library(gplots) #for heatmap.2
onco <- read.csv("oncosig.up.csv", sep=",", header=T, stringsAsFactors = F)
oncolist <- colnames(onco)

pin <- read.csv("../human.pin.csv", header=T, stringsAsFactors = F)
geneA <- pin$geneA
geneB <- pin$geneB

cancer.dim <- length(oncolist)
cancer.matrix <- matrix(0, nrow=cancer.dim, ncol=cancer.dim)
for (i in 1:cancer.dim) {
  for (j in 1:cancer.dim) {
    A.genes <- onco[[i]]
    B.genes <- onco[[j]]
    #number of interactions between cancerA and B
    AB.int <- length(which(((geneA %in% A.genes) & (geneB %in% B.genes)) |
                             (geneA %in% B.genes) & (geneB %in% A.genes)))
    cancer.matrix[i, j] = cancer.matrix[i, j] + AB.int
  }
}

write.table(cancer.matrix, file="oncosig.int.csv", sep=",", col.names=F,
            row.names=F, quote=F)

colnames(cancer.matrix) <- oncolist
rownames(cancer.matrix) <- oncolist

my_palette <- colorRampPalette(c("blue2", "white", "red2"))(n = 29)
colors= c(seq(0,99, length=10), seq(100, 199, length=10), seq(200,300, length=10))


png("onco-onco-int.heatmap.png", width=10, height=9, res=600, units="in")
heatmap.2(cancer.matrix, col=my_palette, trace='none', breaks=colors, 
          key.xlab=NA, key.title="Interaction Count", key.ylab=NA, key.xtickfun = NULL, key.ytickfun = NULL,
          keysize=1,
          srtCol=45, adjCol=c(1,0), dendrogram = "both",
          margins=c(9,10), sepwidth=c(0.01,0.01), 
          sepcolor="grey", colsep=1:cancer.dim, rowsep=1:cancer.dim)
dev.off()

