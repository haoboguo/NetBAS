# overlaps of gene sets
library('gplots')
list.file <- read.csv("list", header=F, stringsAsFactors=F)
hmk.set.list <- list.file$V1

dims <- length(hmk.set.list)

om <- matrix(0, ncol = dims, nrow = dims)

for (i in 1:dims) {
  ith.name <- paste(hmk.set.list[i], "/", hmk.set.list[i], ".csv", sep="")
  ith.file <- read.csv(ith.name, header=T, stringsAsFactors=F)
  ith.set <- ith.file$gene
  for (j in 1:dims) {
    jth.name <- paste(hmk.set.list[j], "/", hmk.set.list[j], ".csv", sep="")
    jth.file <- read.csv(jth.name, header = T, stringsAsFactors = F)
    jth.set <- jth.file$gene
    overlap <- length(which(jth.set %in% ith.set))
    om[i,j] <- om[i,j] + overlap
  }
}

colnames(om) <- hmk.set.list
rownames(om) <- hmk.set.list

write.table(om, file = "Overlap.Matrix.tab", row.names = T, col.names = T, quote = T)

colors = c(seq(0,15,length=10),seq(16,30,length=10),seq(31,200,length=10))
my_palette <- colorRampPalette(c("white", "red2"))(n = 29)3
#png(filename = "Hallmark.Sets.Overlap.Genes.png",width=28, height=28, res=1200, unit="in")
pdf("hallmark.overlap.pdf", width=30,height=30, paper='special')
heatmap.2(om, trace='none', cellnote = om, dendrogram='none',
          colsep = 1:50, rowsep = 1:50, sepcolor="lightgrey", sepwidth = c(0.01,0.01),
          breaks = colors, col=my_palette, Rowv=F, Colv = F,
          ylab="Hallmark Sets", xlab="Hallmark Sets",
          adjCol=c(0,0.2), adjRow=c(0,0.2), srtRow=45, srtCol=-45,
          scale="none", symbreaks=F, symm=F, symkey = F,
          margins = c(15,15), key.title =NA, key.xlab=NA, key.ylab=NA,
          cexRow = 0.8, cexCol=0.8)
dev.off()
