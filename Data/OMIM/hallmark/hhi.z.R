#
library(gplots) #for heatmap.2

hhi.z.dat <- read.csv("hhi.z.csv", header=F, stringsAsFactors = F)
hhi.z <- matrix(unlist(hhi.z.dat), nrow=50, ncol=50)

hhi.p.dat <- read.csv("hhi.p.csv", header=F, stringsAsFactors = F)
hhi.p <- matrix(unlist(hhi.p.dat), nrow=50, ncol=50)

hhi.logp <- -log10(hhi.p+0.00001)

colnames(hhi.z) <- c(1:50)
row.names(hhi.z) <- c(1:50)
colnames(hhi.p) <- c(1:50)
row.names(hhi.p) <- c(1:50)

my_palette <- colorRampPalette(c("blue2", "white", "red2"))(n = 20)
colors = c(seq(-10,10,length=21))

png("hhi.z.heatmap.png", width=12, height=11, res=1200, units="in")
heatmap.2(hhi.z, col=my_palette, trace='none', breaks=colors, 
          Rowv = F, Colv = F, revC = T,
          key.xlab=NA, key.title="Association Z-score", key.ylab=NA, key.xtickfun = NULL, key.ytickfun = NULL,
          srtCol=90, adjCol=c(1,0), dendrogram = "none",
          margins=c(14,18.5), sepwidth=c(0.01,0.01), symbreaks = TRUE,
          sepcolor="grey", colsep=1:hallmark.dim, rowsep=1:hallmark.dim)
dev.off()

my_palette <- colorRampPalette(c("blue2", "white", "red2"))(n = 20)
colors = c(seq(-10,10,length=21))

png("hhi.z.heatmap.colorbar.png", width=12, height=4, res=1200, units="in")
heatmap.2(hhi.z, col=my_palette, trace='none', breaks=colors, 
          Rowv = F, Colv = F, revC = T,
          key.xlab=NA, key.title="Association Z-score", key.ylab=NA, key.xtickfun = NULL, key.ytickfun = NULL,
          srtCol=90, adjCol=c(1,0), dendrogram = "none",
          margins=c(14,18.5), sepwidth=c(0.01,0.01), symbreaks = TRUE,
          sepcolor="grey", colsep=1:hallmark.dim, rowsep=1:hallmark.dim)
dev.off()

my_palette <- colorRampPalette(c("white", "red2"))(n = 10)
colors = c(seq(0,5,length=11))
png("hhi.p.heatmap.png", width=12, height=11, res=1200, units="in")
heatmap.2(hhi.logp, col=my_palette, trace='none', breaks=colors, 
          Rowv = F, Colv = F, revC = T,
          key.xlab=NA, key.title="Z-score", key.ylab=NA, key.xtickfun = NULL, key.ytickfun = NULL,
          srtCol=90, adjCol=c(1,0), dendrogram = "none",
          margins=c(14,18.5), sepwidth=c(0.01,0.01), symbreaks = F,
          sepcolor="grey", colsep=1:hallmark.dim, rowsep=1:hallmark.dim)
dev.off()

gene.count <- read.csv("gene.count.csv", header=T, stringsAsFactors = F)
v.size <- gene.count$gene.number
p.net <- graph.adjacency(hhi.logp, mode="undirected", weighted = T, diag = F)
#apparently every two sets are connected, we only plot those with more than 600 connections
E(p.net)$weight <- ifelse(abs(E(p.net)$weight) > 4.99, E(p.net)$weight, -1)
pdf("hhi.p.network.pdf", height=8, width=8, paper='special')
plot.igraph(p.net, vertex.label=V(p.net)$name, layout=layout_in_circle,
            edge.color="red", edge.width=E(p.net)$weight/4,vertex.size=sqrt(v.size)*1.2)
dev.off()

my_palette <- colorRampPalette(c("white", "red2"))(n = 10)
colors = c(seq(0,5,length=11))
png("hhi.p.colorbar.png", width=12, height=4, res=1200, units="in")
heatmap.2(hhi.logp, col=my_palette, trace='none', breaks=colors, 
          Rowv = F, Colv = F, revC = T,
          key.xlab=NA, key.title=expression(-log[10](p-value)),
          key.ylab=NA, key.xtickfun = NULL, key.ytickfun = NULL,
          srtCol=90, adjCol=c(1,0), dendrogram = "none",
          margins=c(14,18.5), sepwidth=c(0.01,0.01), symbreaks = F,
          sepcolor="grey", colsep=1:hallmark.dim, rowsep=1:hallmark.dim)
dev.off()

list.file <- read.csv("../../MSigDB.go.pathway/list", header=F, stringsAsFactors = F)
hallmark.name <- list.file$V1
overlap.matrix <- matrix(0, ncol=50, nrow=50)
for (i in 1:50) {
  file.name <- paste("../../MSigDB.go.pathway/", hallmark.name[i], "/",
                     hallmark.name[i], ".csv", sep="")
  file.a <- read.csv(file.name, header=T)
  list.a <- file.a$gene
  for (j in 1:50) {
    file.nameb <- paste("../../MSigDB.go.pathway/", hallmark.name[j], "/",
                        hallmark.name[j], ".csv", sep="")
    file.b <- read.csv(file.nameb, header=T)
    list.b <- file.b$gene
    jac <- length(which(list.a %in% list.b))/length(unique(c(list.a, list.b)))
    overlap.matrix[i,j] = overlap.matrix[i,j] + jac
  }
}

colnames(overlap.matrix) <- c(1:50)
row.names(overlap.matrix) <- c(1:50)
my_palette <- colorRampPalette(c("white", "red2"))(n = 10)
colors = c(seq(0,0.1,length=11))
png("hallmark.jaccard.png", width=12, height=11, res=1200, units="in")
heatmap.2(overlap.matrix, col=my_palette, trace='none', breaks=colors, 
          Rowv=F, Colv=F, revC=T,
          key.xlab=NA, key.title="Jaccard Index", 
          key.ylab=NA, key.xtickfun = NULL, key.ytickfun = NULL,
          srtCol=90, adjCol=c(1,0), dendrogram = "none",
          margins=c(14.5,18), sepwidth=c(0.01,0.01), #symbreaks = TRUE,
          sepcolor="grey", colsep=1:50, rowsep=1:50)
dev.off()

  png("hallmark.jaccard.colorbar.png", width=12, height=4, res=1200, units="in")
  heatmap.2(overlap.matrix, col=my_palette, trace='none', breaks=colors, 
            Rowv=F, Colv=F, revC=T,
            key.xlab=NA, key.title="Jaccard Index", 
            key.ylab=NA, key.xtickfun = NULL, key.ytickfun = NULL,
            srtCol=90, adjCol=c(1,0), dendrogram = "none",
            margins=c(14.5,18), sepwidth=c(0.01,0.01), #symbreaks = TRUE,
            sepcolor="grey", colsep=1:50, rowsep=1:50)
  dev.off()

cut <- read.csv("cutoff/cancer.z.csv", header=F, stringsAsFactors = F)
pdf("cutoff/cancer.reshuffle.pdf", width=6, height=5, paper='special')
hist(cut$V1[2:1001], main="", xlab="", ylab="", breaks=10,
     col="lightgrey", xlim=c(-2,8), ylim=c(0,240))
title(ylab="Frequency", line=2.5, cex.lab=1.2)
title(xlab="Z-score", line=2.5, cex.lab=1.2)
arrows(cut$V1[1], 230, cut$V1[1], 0, length=0.1, col="red2", lwd=2)
text(cut$V1[1], 240, "Original Pair", col="red2")
text(3, 240, "Reshuffled Pairs", col="black")
dev.off()

rdmcut <- read.csv("cutoff/rdm.z.csv", header=F, stringsAsFactors = F)
pdf("cutoff/random.reshuffle.pdf", width=6, height=5, paper='special')
hist(rdmcut$V1[2:1001], main="", xlab="", ylab="", breaks=10,
     col="lightgrey", xlim=c(-2.5,2.5), ylim=c(0,300))
title(ylab="Frequency", line=2.5, cex.lab=1.2)
title(xlab="Z-score", line=2.5, cex.lab=1.2)
arrows(0.5, 265, rdmcut$V1[1], 0, length=0.1, col="red2", lwd=2)
text(0.65, 275, "Original Pair", col="red2")
text(-1.25, 275, "Reshuffled Pairs", col="black")
dev.off()
