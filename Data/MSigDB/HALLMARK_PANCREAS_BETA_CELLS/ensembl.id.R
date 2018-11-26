##
# loading the libraries
library(biomaRt)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(GenomicFeatures)

gene.file <- read.csv("HALLMARK_PANCREAS_BETA_CELLS.csv", header=T, stringsAsFactors=F)
genes <- gene.file$gene

ensembl.id <- genes(EnsDb.Hsapiens.v86,
                    filter=list(GeneNameFilter(genes),GeneIdFilter("ENSG", "startsWith")),
                    return.type="data.frame", columns=c("gene_id"))
write.table(ensembl.id, file="HALLMARK_PANCREAS_BETA_CELLS.geneset.csv", sep=",", row.names = F, col.names = F, quote = T)
