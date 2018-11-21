####
#library('data.table')
rm(list=ls())
human.pathway <- read.csv("human.pathways.table.txt", sep="\t", header=T, stringsAsFactors = F)
pathway <- human.pathway$pathway
resource <- human.pathway$source
gene <- human.pathway$gene

pair <- c("Source", "Pathway", "Gene")

for (i in 1:length(pathway)) {
  path <- as.character(pathway[i])
  reso <- as.character(resource[i])
  gene.expand <- as.character(unlist(strsplit(as.character(gene[i]), ",")))
     for (j in 1:length(gene.expand)) {
        pair <- rbind(pair, c(reso, path, gene.expand[j]))
  }
}

write.table(pair, file="human.pathway.csv", sep="\t", row.names = F, col.names=F, quote=T)
