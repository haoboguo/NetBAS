rm(list=ls())
library(igraph)

human.pin <- read.csv("../../human.pin.csv", header=T, stringsAsFactors = F)
geneA <- human.pin$geneA
geneB <- human.pin$geneB

#human.pin <- read.csv("../../../ms02star/human/ms02.1.csv", header=T, stringsAsFactors=F)
#geneA <- human.pin$id1
#geneB <- human.pin$id2

#All genes appear in PIN
all.list <- unique(c(geneA, geneB))

#the Gene list
geneset <- read.csv("HALLMARK_PI3K_AKT_MTOR_SIGNALING.csv", header=T, stringsAsFactors = F)
gene.list <- geneset$gene

go.file <- read.csv("../../human.all.go.csv", header=T, stringsAsFactors = F)
go.id <- go.file$goid
go.gene <- go.file$gene

enriched <- read.csv("HMK.bp.new.txt", sep="\t", header=T, stringsAsFactors=F)
term <- enriched$GO.ID[2]

term.list <- go.gene[which(go.id %in% term)]

neiborA <- geneA[which(geneB %in% gene.list)]
subA <- geneB[which(geneB %in% gene.list)]
neiborB <- geneB[which(geneA %in% gene.list)]
subB <- geneA[which(geneA %in% gene.list)]

selA <- neiborA[which(neiborA %in% term.list)]
selA2 <- subA[which(neiborA %in% term.list)]

selB <- subB[which(neiborB %in% term.list)]
selB2 <- neiborB[which(neiborB %in% term.list)]

web1 <- cbind(selA, selA2)
web2 <- cbind(selB, selB2)
web <- rbind(web1,web2)

neibor.web <- data.frame(web)
neibor.graph <- graph.data.frame(neibor.web, directed = F)

'%ni%' <- Negate('%in%')
blue.id <- which(as_ids(V(neibor.graph)) %in% term.list & as_ids(V(neibor.graph)) %ni% gene.list)
red.id <- which(as_ids(V(neibor.graph)) %in% gene.list & as_ids(V(neibor.graph)) %ni% term.list)
green.id <- which(as_ids(V(neibor.graph)) %in% gene.list & as_ids(V(neibor.graph)) %in% term.list)

color <- rep("NA", times=length(V(neibor.graph)))
color[red.id] <- rep("red", times=length(red.id))
color[blue.id] <- rep("blue", times=length(blue.id))
color[green.id] <- rep("green", times=length(green.id))
V(neibor.graph)$color <- color

pdf("HALLMARK_PI3K_AKT_MTOR_SIGNALING.top2.neibor.pdf")
plot.igraph(neibor.graph, vertex.color=V(neibor.graph)$color, vertex.size=5, vertex.label=NA, edge.width=1)
dev.off()

#neibor.cluster <- cluster_walktrap(neibor.graph)

#plot(neibor.cluster, neibor.graph, modularity=T, vertex.label=NA, vertex.size=5)
