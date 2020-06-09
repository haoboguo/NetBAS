#
library(gplots) #for heatmap.2

list.file <- read.csv("../../MSigDB.go.pathway/list", header=F, stringsAsFactors = F)
list <- list.file$V1
hallmark.dim <- 50

hhi.dat <- read.csv("hhi.z.csv", header=F, stringsAsFactors = F)
hhi.z <- matrix(unlist(hhi.dat), nrow=hallmark.dim, ncol=hallmark.dim)

colnames(hhi.z) <- list
rownames(hhi.z) <- list

my_palette <- colorRampPalette(c("blue2", "white", "red2"))(n = 20)
#colors= c(seq(-5,-1.5, length=10), seq(-1.49, 1.49, length=10), seq(1.5,5, length=10))
colors = c(seq(-10,10,length=21))

png("hhi.z.heatmap.png", width=12, height=11, res=600, units="in")
heatmap.2(hhi.z, col=my_palette, trace='none', breaks=colors, 
          key.xlab=NA, key.title="Z-score", key.ylab=NA, key.xtickfun = NULL, key.ytickfun = NULL,
          srtCol=45, adjCol=c(1,0.5), dendrogram = "both",
          margins=c(14,18.5), sepwidth=c(0.01,0.01), symbreaks = TRUE,
          sepcolor="grey", colsep=1:hallmark.dim, rowsep=1:hallmark.dim)
dev.off()

hhi.z.mat <- as.matrix(hhi.z)
hhi.z.net <- graph.adjacency(hhi.z.mat, mode="undirected", weighted=T, diag=F)

top2 <- quantile(hhi.z, probs=seq(0,1,1/50))[50]
bottom2 <- quantile(hhi.z, probs=seq(0,1,1/50))[2]

E(hhi.z.net)$color <- ifelse(E(hhi.z.net)$weight > 0, "red", "blue")
coloring <- E(hhi.z.net)$color
hhi.weight <- ifelse((E(hhi.z.net)$weight > top2 | (E(hhi.z.net)$weight < bottom2)), 
                     abs(E(hhi.z.net)$weight), -0.5)

pdf("hhi.z.network.pdf", height=8, width=8, paper='special')
plot.igraph(hhi.z.net, vertex.label=c(1:50), layout=layout_in_circle,
            edge.color = coloring, edge.width=hhi.weight)
dev.off()
