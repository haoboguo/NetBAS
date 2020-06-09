discurate <- read.csv("curated_gene_disease_associations.tsv", 
                      header=T, stringsAsFactors = F, sep="\t")

id.list <- unique(discurate$diseaseId[which(discurate$diseaseType == "disease")])

name.list <- c()
for (k in 1:length(id.list)) {
  name.tmp <- unique(discurate$diseaseName[which(discurate$diseaseId == id.list[k])])
  name.list <- rbind(name.list, name.tmp)
}

out.data <- c("gene", "ID", "Name")

for (i in 1:length(id.list)) {
  scores <- discurate$score[which(discurate$diseaseId == id.list[i])]
  if (length(scores) > 200) {
    scores <- sort(scores)
    cutoff <- scores[length(scores)-200]
    gene.list <- discurate$geneSymbol[which((discurate$diseaseId == id.list[i]) &
                                        (discurate$score > cutoff))]
    for (j in 1:length(gene.list)) {
      out.data <- rbind(out.data, c(gene.list[j], id.list[i], name.list[i]))
    }
  } else {
  gene.list <- discurate$geneSymbol[which(discurate$diseaseId == id.list[i])]
  for (j in 1:length(gene.list)) {
    out.data <- rbind(out.data, c(gene.list[j], id.list[i], name.list[i]))
  }
  }
}

write.table(out.data, file="Discurate.full.csv", sep=",", col.names=F, row.names = F, quote=T)