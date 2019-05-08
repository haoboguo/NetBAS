rm(list=ls())
library(igraph)

human.pin <- read.csv("data/human.pin.csv", header=T, stringsAsFactors = F)
geneA <- human.pin$geneA
geneB <- human.pin$geneB

#All genes appear in PIN
all.list <- unique(c(geneA, geneB))

#the Gene list
geneset <- read.csv("genes.csv", header=T, stringsAsFactors = F)
gene.list <- geneset$gene

go.file <- read.csv("data/human.go.csv", header=T, stringsAsFactors = F)
go.id <- go.file$goid
go.gene <- go.file$gene

go.type <- c("bp", "cc", "mf")

#loop1, different GO types
for (i in 1:3) {

gofile.name <- paste("output/", "genes.", go.type[i], ".txt", sep="")

enriched <- read.csv(gofile.name, sep="\t", header=T, stringsAsFactors=F)

#loop2, top 5 enriched terms
for (k in 1:5) {
term <- enriched$GO.ID[k]
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
yellow.id <- which(as_ids(V(neibor.graph)) %in% gene.list & as_ids(V(neibor.graph)) %in% term.list)

color <- rep("NA", times=length(V(neibor.graph)))
color[red.id] <- rep("red", times=length(red.id))
color[blue.id] <- rep("blue", times=length(blue.id))
color[yellow.id] <- rep("yellow", times=length(yellow.id))
V(neibor.graph)$color <- color

pdf.name <- paste("network/", "genes.", go.type[i], ".", k, ".neibor.pdf", sep="")
pdf(pdf.name)
plot.igraph(neibor.graph, vertex.color=V(neibor.graph)$color, vertex.size=10, vertex.label=NA, edge.width=1)
dev.off()
}#loop1
}#loop2

