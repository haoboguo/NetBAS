sel <- read.csv("DisGeNET.selected.csv", header=T, stringsAsFactors = F)
list <- unique(sel$ID)

set.size <- c("ID", "Name", "Size")
for (k in 1:length(list)) {
  size <- length(sel$gene[which(sel$ID == list[k])])
  name <- unique(sel$Name[which(sel$ID == list[k])])
  set.size <- rbind(set.size, c(list[k], name, size))
}

pin <- read.csv("../human.pin.csv", header=T, stringsAsFactors = F)
geneA <- pin$geneA
geneB <- pin$geneB

dim <- length(list)

disorder.matrix <- matrix(0, nrow=dim, ncol=dim)

for (i in 1:dim) {
  for (j in 1:dim) {
    disA <- list[i]
    disB <- list[j]
    A.genes <- sel$gene[which(sel$ID == disA)]
    B.genes <- sel$gene[which(sel$ID == disB)]
    #number of interactions between cancerA and B
    AB.int <- length(which(((geneA %in% A.genes) & (geneB %in% B.genes)) |
                             (geneA %in% B.genes) & (geneB %in% A.genes)))
    disorder.matrix[i, j] = disorder.matrix[i, j] + AB.int
  }
}

write.table(disorder.matrix, file="ddi.csv", sep=",", col.names=F,
            row.names=F, quote=F)

colnames(disorder.matrix) <- list
rownames(disorder.matrix) <- list

my_palette <- colorRampPalette(c("blue2", "white", "red2"))(n = 19)
colors= c(seq(0, 4, length=20))


png("ddi.log10.heatmap.png", width=10, height=9, res=600, units="in")
heatmap.2(log10(disorder.matrix+1), col=my_palette, trace='none', breaks=colors, 
          key.xlab=NA, key.title="Interaction Count", key.ylab=NA, key.xtickfun = NULL, key.ytickfun = NULL,
          keysize=1,
          srtCol=45, adjCol=c(1,0.5), dendrogram = "both",
          margins=c(8,12), sepwidth=c(0.01,0.01), 
          sepcolor="grey", colsep=1:dim, rowsep=1:dim)
dev.off()