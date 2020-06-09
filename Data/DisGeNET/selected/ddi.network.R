#
library(gplots) #for heatmap.2
library(igraph)

sel <- read.csv("DisGeNET.selected.csv", header=T, stringsAsFactors = F)
list <- unique(sel$ID)
#list <- c(1:50)
disorder.dim <- 50

ddi.dat <- read.csv("ddi.z.csv", header=F, stringsAsFactors = F)
ddi.z <- matrix(unlist(ddi.dat), nrow=disorder.dim, ncol=disorder.dim)

colnames(ddi.z) <- list
rownames(ddi.z) <- list

my_palette <- colorRampPalette(c("blue2", "white", "red2"))(n = 20)
#colors= c(seq(-5,-1.5, length=10), seq(-1.49, 1.49, length=10), seq(1.5,5, length=10))
colors = c(seq(-10,10,length=21))

png("ddi.z.heatmap.png", width=12, height=11, res=600, units="in")
heatmap.2(ddi.z, col=my_palette, trace='none', breaks=colors, 
          key.xlab=NA, key.title="Z-score", key.ylab=NA, key.xtickfun = NULL, key.ytickfun = NULL,
          #srtCol=45, adjCol=c(1,0), dendrogram = "both",
          margins=c(14,18.5), sepwidth=c(0.01,0.01), symbreaks = TRUE,
          sepcolor="grey", colsep=1:disorder.dim, rowsep=1:disorder.dim)
dev.off()

ddi.z.mat <- as.matrix(ddi.z)
ddi.z.net <- graph.adjacency(ddi.z.mat, mode="undirected", weighted=T, diag=F)

top2 <- quantile(ddi.z, probs=seq(0,1,1/10))[10]
bottom2 <- quantile(ddi.z, probs=seq(0,1,1/10))[2]

E(ddi.z.net)$color <- ifelse(E(ddi.z.net)$weight > 0, "red", "blue")
coloring <- E(ddi.z.net)$color
ddi.weight <- ifelse((E(ddi.z.net)$weight > top2 | (E(ddi.z.net)$weight < bottom2)), 
                     abs(E(ddi.z.net)$weight), -0.5)
#ddi.weight <- abs(E(ddi.z.net)$weight)

pdf("ddi.z.network.pdf", height=8, width=8, paper='special')
plot.igraph(ddi.z.net, vertex.label=c(1:50), #layout=layout_in_circle,
            edge.color = coloring, edge.width=ddi.weight, vertex.size=10)
dev.off()

