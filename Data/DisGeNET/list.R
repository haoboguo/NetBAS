discurate <- read.csv("Discurate.full.non.MIR.csv", 
                      header=T, stringsAsFactors = F, sep=",")

id.list <- unique(discurate$ID)

name.list <- c()
for (k in 1:length(id.list)) {
  name.tmp <- unique(discurate$Name[which(discurate$ID == id.list[k])])
  name.list <- rbind(name.list, name.tmp)
}

pin <- read.csv("../human.pin.csv", header=T, stringsAsFactors = F)
geneA <- pin$geneA
geneB <- pin$geneB

pin.list <- unique(c(geneA, geneB))

list.data <- c("ID", "Name", "GN")
for (j in 1:length(id.list)) {
  disease.id <- id.list[j]
  disease.name <- name.list[j]
  disease.gene <- discurate$gene[which(discurate$ID == disease.id)]
  overlap <- length(which(disease.gene %in% pin.list))
  list.data <- rbind(list.data, c(disease.id, disease.name, overlap))
}

head(list.data)

write.table(list.data, file="all.ids.csv", col.names = F, row.names = F, quote = T, sep=",")

test <- read.csv("all.ids.csv", header=T, stringsAsFactors = F)
length(which(test$GN > 14))
