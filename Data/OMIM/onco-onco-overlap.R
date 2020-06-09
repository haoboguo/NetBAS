onco <- read.csv("onco.signature.uponly.csv", sep=",", header=T, stringsAsFactors = F)
oncolist <- colnames(onco)
onco[onco == ""] <- NA
cancer.dim <- length(oncolist)

overlap.dim <- length(oncolist)
overlap.matrix <- matrix(0, nrow=cancer.dim, ncol=cancer.dim)
for (i in 1:overlap.dim) {
  for (j in 1:overlap.dim) {
    A.genes <- onco[[i]][!is.na(onco[[i]])]
    B.genes <- onco[[j]][!is.na(onco[[j]])]
    #number of overlaps between cancerA and B
    AB.over <- length(which(A.genes %in% B.genes))
    overlap.matrix[i, j] = overlap.matrix[i, j] + AB.over
  }
}

write.table(overlap.matrix, file="onco-onco.overlap.csv", sep=",", col.names=F,
            row.names=F, quote=F)

colnames(overlap.matrix) <- oncolist
rownames(overlap.matrix) <- oncolist

my_palette <- colorRampPalette(c("blue2", "white", "red2"))(n = 9)
colors= c(seq(0,18, length=10))

png("onco-onco-overlap.heatmap.png", width=12, height=11, res=600, units="in")
heatmap.2(overlap.matrix, col=my_palette, trace='none', breaks=colors, 
          key.xlab=NA, key.title="Overlap Gene Numbers", key.ylab=NA, key.xtickfun = NULL, key.ytickfun = NULL,
          keysize=1,
          srtCol=45, adjCol=c(1,0), dendrogram = "both",
          margins=c(9,12), sepwidth=c(0.01,0.01), 
          sepcolor="grey", colsep=1:cancer.dim, rowsep=1:cancer.dim)
dev.off()

#hm <- heatmap.2(overlap.matrix, dendrogram = "both")
#hc.row <- as.hclust(hm$rowDendrogram)

#pdf("onco-onco-overlap-dendrogram.pdf", width=10, height=6,paper='special')
#plot(hc.row, main="Oncology Signature Genes Dendrogram", cex=.8)
#dev.off()
