rm(list=ls())
library(WebGestaltR)
refFile <- system.file("extdata", "referenceGenes.txt", package="WebGestaltR")
WebGestaltR(enrichMethod="ORA", organism="hsapiens",
    enrichDatabase="pathway_KEGG", interestGeneFile="HALLMARK_XENOBIOTIC_METABOLISM.geneset.txt",
    interestGeneType="ensembl_gene_id", referenceGeneFile=refFile,
    referenceGeneType="genesymbol", isOutput=TRUE,
    outputDirectory=getwd(), projectName=NULL)
