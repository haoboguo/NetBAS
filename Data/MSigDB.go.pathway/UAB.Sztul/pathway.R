rm(list=ls())
library(igraph)

human.pin <- read.csv("../../human.pin.csv", header=T, stringsAsFactors = F)
geneA <- human.pin$geneA
geneB <- human.pin$geneB

#the Gene list
geneset <- read.csv("genes.csv", header=T, stringsAsFactors = F)
gene.list <- geneset$gene

pw.file <- read.csv("../../human.pathway.ids.csv", sep="\t", header=T, stringsAsFactors=F)
pw.gene <- pw.file$Gene
pw.id <- pw.file$ID
pw.pathway <- pw.file$Pathway
pw.source <- pw.file$Source

pw.type <- c("kegg", "biocarta", "reactome")

#loop1
for (i in 1:3) {

pwfile.name <- paste("HMK.", pw.type[i], ".new.txt", sep="")

enriched <- read.csv(pwfile.name, sep="\t", header=T, stringsAsFactors=F)

#loop2
for (k in 1:10) {
term <- enriched$ID[k]
term.list <- pw.gene[which(pw.id %in% term)]

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

pdf.name <- paste("UAB.Sztul.", pw.type[i], ".", k, ".neibor.pdf", sep="")
pdf(pdf.name)
plot.igraph(neibor.graph, vertex.color=V(neibor.graph)$color, vertex.size=10, vertex.label=NA, edge.width=1)
dev.off()
}#loop1
}#loop2

