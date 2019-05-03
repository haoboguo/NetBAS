rm(list=ls())

#human.pin <- read.csv("../../human.pin.csv", header=T, stringsAsFactors = F)
#geneA <- human.pin$geneA
#geneB <- human.pin$geneB

human.pin <- read.csv("../../../ms02star/human/ms02.1.csv", header=T, stringsAsFactors=F)
geneA <- human.pin$id1
geneB <- human.pin$id2

#All genes appear in PIN
#all.list <- unique(c(geneA, geneB))

#the Gene list
geneset <- read.csv("HALLMARK_HEDGEHOG_SIGNALING.csv", header=T, stringsAsFactors = F)
gene.list <- geneset$gene

go.file <- read.csv("../../human.all.go.csv", header=T, stringsAsFactors = F)
go.id <- go.file$goid
go.gene <- go.file$gene

term <- c("GO:1902287")

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

#V(neibor.graph)$color <- "blue"
#V(neibor.graph)$color[1:4] <- "red"

V(neibor.graph)$color <- c("green", "blue", "blue", "red")
pdf("go.1902287.ms02.pdf", width=3, height=3, paper='special')
plot.igraph(neibor.graph, vertex.color=V(neibor.graph)$color, vertex.size=30, edge.width=2)#, vertex.label=NA)
dev.off()

neibor.cluster <- cluster_walktrap(neibor.graph)

plot(neibor.cluster, neibor.graph, modularity=T, vertex.label=NA, vertex.size=5)
