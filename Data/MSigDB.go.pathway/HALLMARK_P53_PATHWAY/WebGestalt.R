rm(list=ls())
library(WebGestaltR)
WebGestaltR(enrichMethod="NTA", organism="hsapiens",
    enrichDatabase="network_PPI_BIOGRID", interestGeneFile="HALLMARK_P53_PATHWAY.geneset.txt",
    interestGeneType="ensembl_gene_id",sigMethod="fdr", 
    outputDirectory=getwd(), highlightType = 'Neighbors', neighborNum = 10,
    networkConstructionMethod="Network_Expansion")
