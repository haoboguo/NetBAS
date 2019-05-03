detach("package:pathfindR", unload=TRUE)
detach("package:pathview", unload=TRUE)
library(pathfindR)
library(stringr)

gene.file <- read.csv("HALLMARK_REACTIVE_OXIGEN_SPECIES_PATHWAY.csv", header=T, stringsAsFactors=F)
gene.list <- gene.file$gene
input <- data.frame(gene.list)
input$change <- as.numeric(rep(c(0), times=length(gene.list)))
input$pVal <- as.numeric(rep(c(0.05), times=length(gene.list)))

run_pathfindR(input, p_val_threshold=0.05, visualize_pathways=T,
              human_genes=T, enrichment_threshold=0.05,
              adj_method="bonferroni", search_method="GR")
