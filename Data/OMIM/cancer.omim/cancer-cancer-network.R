library(igraph)

cci <- read.csv("cci.csv", header=F, stringsAsFactors = F)
cancer <- read.csv("cancer.gene.csv", header=T, stringsAsFactors = F)
cancer.name <- unique(cancer$NAME[!is.na(cancer$SN)])
cancer.dim <- length(cancer.name)

colnames(cci) <- cancer.name
row.names(cci) <- cancer.name

cci.mat <- as.matrix(cci)

cci.net <- graph.adjacency(cci.mat, mode="undirected", weighted = T, diag = F)
summary(cci.net)

E(cci.net)$weight

pdf("cancer-cancer.network.pdf", height=8, width = 8, paper='special')
plot.igraph(cci.net, vertext.label=V(cci.net)$name, layout=layout_in_circle,
            edge.color="lightblue", edge.width=E(cci.net)$weight/2)
dev.off()

overlap <- read.csv("cancer-cancer.overlap.csv", header=F, stringsAsFactors = F)
colnames(overlap) <- cancer.name
row.names(overlap) <- cancer.name

overlap.mat <- as.matrix(overlap)
overlap.net <- graph.adjacency(overlap.mat, mode="undirected", weighted=T, diag=F)

plot.igraph(overlap.net, vertext.label=V(overlap.net)$name, layout=layout_in_circle,
            edge.color="lightblue", edge.width=E(overlap.net)$weight)

cci.z <- read.csv("cci.z.csv", header=F, stringsAsFactors = F)
colnames(cci.z) <- cancer.name
row.names(cci.z) <- cancer.name
cci.z.mat <- as.matrix(cci.z)
cci.z.net <- graph.adjacency(cci.z.mat, mode="undirected", weighted=T, diag=F)
summary(cci.z.net)

E(cci.z.net)$color <- ifelse(E(cci.z.net)$weight > 0, "red", "blue")
coloring <- E(cci.z.net)$color
E(cci.z.net)$weight <- ifelse(E(cci.z.net)$weight > 0, E(cci.z.net)$weight, -E(cci.z.net)$weight)

pdf("cci.z.network.pdf", height=8, width=8, paper='special')
plot.igraph(cci.z.net, vertext.label=V(cci.z.net)$name, layout=layout_in_circle,
            edge.color = coloring, edge.width=E(cci.z.net)$weight)
dev.off()