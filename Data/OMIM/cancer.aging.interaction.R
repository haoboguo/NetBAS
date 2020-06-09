#
library(gplots) #for heatmap.2
cancer <- read.csv("cancer.gene.csv", header=T, stringsAsFactors = F)
cancer.name <- unique(cancer$NAME[!is.na(cancer$SN)])

aging <- read.csv("microarray.csv", header=T, stringsAsFactors = F)

over.list <- aging$Gene[which(aging$Expression == "Overexpressed")]
under.list <- aging$Gene[which(aging$Expression == "Underexpressed")]

pin <- read.csv("../human.pin.csv", header=T, stringsAsFactors = F)

geneA <- pin$geneA
geneB <- pin$geneB

#calculate association between a pair of two cancers

cancer.dim <- length(cancer.name)

cancer.aging.matrix <- matrix(0, ncol=cancer.dim, nrow=2)

for (i in 1:cancer.dim) {
    cancer.gene <- cancer.name[i]
    cancer.list <- cancer$Gene[which(cancer$NAME == cancer.gene)]
    #number of interactions between cancerA and over
    CO.int <- length(which(((geneA %in% cancer.list) & (geneB %in% over.list)) |
                            (geneA %in% over.list) & (geneB %in% cancer.list)))
    cancer.aging.matrix[1, i] <- cancer.aging.matrix[1, i] + CO.int
    CU.int <- length(which(((geneA %in% cancer.list) & (geneB %in% under.list)) |
                             (geneA %in% under.list) & (geneB %in% cancer.list)))
    cancer.aging.matrix[2, i] <- cancer.aging.matrix[2, i] + CU.int    
}

write.table(cancer.aging.matrix, file="cancer.aging.interaction.csv", sep=",", col.names=F,
            row.names=F, quote=F)

colnames(cancer.aging.matrix) <- cancer.name
rownames(cancer.aging.matrix) <- c("Over Expression in Aging", "Under Expression in Aging")

my_palette <- colorRampPalette(c("blue2", "white", "red2"))(n = 29)
colors= c(seq(0,3, length=10), seq(3.1, 6, length=10), seq(6.1,10, length=10))
  

png("cc-aging.int.heatmap.png", width=10, height=4, res=600, units="in")
heatmap.2(cancer.aging.matrix, col=my_palette, trace='none', breaks=colors, 
          key.xlab=NA, key.title="Interaction Count", key.ylab=NA, key.xtickfun = NULL, key.ytickfun = NULL,
          keysize=1,
          srtCol=45, adjCol=c(1,0), dendrogram = "both",
          margins=c(8,12), sepwidth=c(0.01,0.01), cexRow = 0.8, cexCol = 0.8,
          sepcolor="grey", colsep=1:cancer.dim, rowsep=1:cancer.dim)
dev.off()

