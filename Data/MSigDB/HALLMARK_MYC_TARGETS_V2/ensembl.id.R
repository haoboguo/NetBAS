##
# loading the libraries
library(biomaRt)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(GenomicFeatures)

gene.file <- read.csv("HALLMARK_MYC_TARGETS_V2.csv", header=T, stringsAsFactors=F)
genes <- gene.file$gene

ensembl.id <- genes(EnsDb.Hsapiens.v86,
                    filter=list(GeneNameFilter(genes),GeneIdFilter("ENSG", "startsWith")),
                    return.type="data.frame", columns=c("gene_id"))
write.table(ensembl.id, file="HALLMARK_MYC_TARGETS_V2.geneset.txt", sep=",", row.names = F, col.names = F, quote = T)
