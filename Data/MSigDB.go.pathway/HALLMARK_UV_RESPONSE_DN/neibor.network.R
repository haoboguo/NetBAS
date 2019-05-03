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
geneset <- read.csv("HALLMARK_UV_RESPONSE_DN.csv", header=T, stringsAsFactors = F)
gene.list <- geneset$gene

go.file <- read.csv("../../human.all.go.csv", header=T, stringsAsFactors = F)
go.id <- go.file$goid
go.gene <- go.file$gene

go.type <- c("bp", "cc", "mf")

record <- c()
#loop1
for (i in 1:3) {

gofile.name <- paste("HMK.", go.type[i], ".new.txt", sep="")

enriched <- read.csv(gofile.name, sep="\t", header=T, stringsAsFactors=F)

#loop2
for (k in 1:10) {

title <-  paste(go.type[i], ".", k, sep="")

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
green.id <- which(as_ids(V(neibor.graph)) %in% gene.list & as_ids(V(neibor.graph)) %in% term.list)

record <- rbind(record, c(title, length(blue.id), length(red.id), length(green.id)))

}#loop1
}#loop2

write.table(record, file="neighbor.network.dist.csv", sep=",", row.names=F, col.names=F, quote=F)
