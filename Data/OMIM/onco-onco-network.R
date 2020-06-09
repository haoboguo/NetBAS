#
library(gplots) #for heatmap.2

onco <- read.csv("onco.signature.uponly.csv", sep=",", header=T, stringsAsFactors = F)
oncolist <- colnames(onco)
cancer.dim <- length(oncolist)

ooi.dat <- read.csv("ooi.z.csv", header=F, stringsAsFactors = F)
ooi.z <- matrix(unlist(ooi.dat), nrow=cancer.dim, ncol=cancer.dim)

colnames(ooi.z) <- oncolist
rownames(ooi.z) <- oncolist

ooi.z.mat <- as.matrix(ooi.z)
ooi.z.net <- graph.adjacency(ooi.z.mat, mode="undirected", weighted=T, diag=F)

E(ooi.z.net)$color <- ifelse(E(ooi.z.net)$weight > 0, "red", "blue")
coloring <- E(ooi.z.net)$color

#E(ooi.z.net)$weight <- -E(ooi.z.net)$weight
#E(ooi.z.net)$weight <- ifelse(E(ooi.z.net)$weight > 0, E(ooi.z.net)$weight, -E(ooi.z.net)$weight)

pdf("onco-onco-pos-z.pdf", height=8, width=8, paper='special')
plot.igraph(ooi.z.net, vertex.label=NA, layout=layout_in_circle, vertex.size=10,
            edge.width=E(ooi.z.net)$weight/5)
dev.off()

E(ooi.z.net)$weight <- -E(ooi.z.net)$weight

pdf("onco-onco-neg-z.pdf", height=8, width=8, paper='special')
plot.igraph(ooi.z.net, vertex.label=NA, layout=layout_in_circle, vertex.size=10,
            edge.width=E(ooi.z.net)$weight/5)
dev.off()