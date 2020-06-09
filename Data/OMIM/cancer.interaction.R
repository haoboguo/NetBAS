#

library(gplots) #for heatmap.2
cancer <- read.csv("cancer.gene.csv", header=T, stringsAsFactors = F)

cancer.name <- unique(cancer$NAME[!is.na(cancer$SN)])

pin <- read.csv("../human.pin.csv", header=T, stringsAsFactors = F)

geneA <- pin$geneA
geneB <- pin$geneB

#calculate association between a pair of two cancers

cancer.dim <- length(cancer.name)

cancer.matrix <- matrix(0, nrow=cancer.dim, ncol=cancer.dim)

for (i in 1:cancer.dim) {
  for (j in 1:cancer.dim) {
    cancerA <- cancer.name[i]
    cancerB <- cancer.name[j]
    A.genes <- cancer$Gene[which(cancer$NAME == cancerA)]
    B.genes <- cancer$Gene[which(cancer$NAME == cancerB)]
    #number of interactions between cancerA and B
    AB.int <- length(which(((geneA %in% A.genes) & (geneB %in% B.genes)) |
                            (geneA %in% B.genes) & (geneB %in% A.genes)))
    cancer.matrix[i, j] = cancer.matrix[i, j] + AB.int
  }
}

write.table(cancer.matrix, file="cci.csv", sep=",", col.names=F,
            row.names=F, quote=F)

colnames(cancer.matrix) <- cancer.name
rownames(cancer.matrix) <- cancer.name

my_palette <- colorRampPalette(c("blue2", "red2"))(n = 9)
colors= c(seq(0,5, length=10))


png("cc-inter.heatmap.png", width=10, height=9, res=600, units="in")
heatmap.2(cancer.matrix, col=my_palette, trace='none', breaks=colors, 
          key.xlab=NA, key.title="Interaction Count", key.ylab=NA, key.xtickfun = NULL, key.ytickfun = NULL,
          keysize=1,
          srtCol=45, adjCol=c(1,0), dendrogram = "both",
          margins=c(8,12), sepwidth=c(0.01,0.01), 
          sepcolor="grey", colsep=1:cancer.dim, rowsep=1:cancer.dim)
dev.off()

