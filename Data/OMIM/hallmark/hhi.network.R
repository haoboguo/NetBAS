library(igraph)

list.file <- read.csv("../../MSigDB.go.pathway/list", header=F, stringsAsFactors = F)
hhi <- read.csv("hhi.csv", header=F, stringsAsFactors = F)
hallmark.name <- list.file$V1
hallmark.dim <- length(hallmark.name)

colnames(hhi) <- hallmark.name
row.names(hhi) <- hallmark.name

hhi.mat <- as.matrix(hhi)

hhi.net <- graph.adjacency(hhi.mat, mode="undirected", weighted = T, diag = F)
summary(hhi.net)

E(hhi.net)$weight

pdf("hallmark-hallmark.network.pdf", height=8, width = 8, paper='special')
plot.igraph(hhi.net, vertex.label=NA, layout=layout_in_circle,
            edge.color="lightblue", edge.width=log(E(hhi.net)$weight))
dev.off()

overlap <- read.csv("hallmark-hallmark.overlap.csv", header=F, stringsAsFactors = F)
colnames(overlap) <- hallmark.name
row.names(overlap) <- hallmark.name

overlap.mat <- as.matrix(overlap)
overlap.net <- graph.adjacency(overlap.mat, mode="undirected", weighted=T, diag=F)

plot.igraph(overlap.net, vertext.label=V(overlap.net)$name, layout=layout_in_circle,
            edge.color="lightblue", edge.width=E(overlap.net)$weight)

hhi.z <- read.csv("hhi.z.csv", header=F, stringsAsFactors = F)
colnames(hhi.z) <- hallmark.name
row.names(hhi.z) <- hallmark.name
hhi.z.mat <- as.matrix(hhi.z)
hhi.z.net <- graph.adjacency(hhi.z.mat, mode="undirected", weighted=T, diag=F)
summary(hhi.z.net)

E(hhi.z.net)$color <- ifelse(E(hhi.z.net)$weight > 0, "red", "blue")
coloring <- E(hhi.z.net)$color
E(hhi.z.net)$weight <- ifelse(abs(E(hhi.z.net)$weight) > 5, abs(E(hhi.z.net)$weight), 0)

pdf("hhi.z.network.5plus.pdf", height=8, width=8, paper='special')
plot.igraph(hhi.z.net, vertex.label=c(1:50), layout=layout_in_circle,
            edge.color = coloring, edge.width=E(hhi.z.net)$weight/4)
dev.off()
