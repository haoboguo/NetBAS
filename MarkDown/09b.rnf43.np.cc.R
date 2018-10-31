## read the original network
library(plyr)
network <- read.csv("../Data/human.pin.csv", header=T, stringsAsFactors=F)
geneA <- network$geneA
geneB <- network$geneB

## the GO terms for biological processes (BP)
GOcategory.file <- read.csv("../Data/human.cc.term.csv",header=TRUE, stringsAsFactors=F)
cc.go.cat <- GOcategory.file$GO.term
cc.dim <- length(cc.go.cat)

# Gene-GO term file
# Note in human PIN the standard gene names have been used
GOterm.file <- read.csv("../Data/human.cc.gene.term.csv", header=T, stringsAsFactors=F)
cc.GO.gene <- GOterm.file$gene  #it should be changed to System for yeast pin
cc.GO.term <- GOterm.file$GO.term

#The gene list is selected from the RNF43 set 
## Steinhart, Z.; et al. Genome-wide CRISPR screens reveal a Wnt-FZD5 signaling circut
## as a druggable vulnerability of RNF43-mutant pancreatic tumors. Nature Medicine 2017.
## BF factors are used here (Table S4), the mean Z-scores of three pancreatic tumor cell lines
## Quantile 20 indicates 858 genes with the top 5% Z-scores (including FZD5, WLS, etc.)
panc.file <- read.csv("../rnf43.csv",header=TRUE,stringsAsFactors=F)
panc.gene <- panc.file$gene
panc.panc <- panc.file$panc.mean
quant <- quantile(panc.panc, probs = seq(0,1,1/20))

all.count <- c()

#gene lists Q1 to Q20
  genelist.1  <- panc.gene[which(panc.panc <= quant[2])]
  glA <- geneB[which(geneA %in% genelist.1)] 
  glB <- geneA[which(geneB %in% genelist.1)]
  gl.all <- c(glA, glB)
  gene.count <- count(gl.all)
  orfs <- gene.count$x
  freqs <- gene.count$freq
  vec <- numeric(length=cc.dim)

#Note that this loop is necessary to avoid miss counting owing to the redundancy of genes
  for (i in 1:length(orfs)) {
    orf.ith <- orfs[i]
    orf.freq <- freqs[i]
    orf.term <- cc.GO.term[which(cc.GO.gene %in% orf.ith)]
    for (j in 1:length(orf.term)) {
      na <- which(cc.go.cat %in% orf.term[j])
      vec[na] <- vec[na] + orf.freq
    } 
  }
  all.count <- c(all.count, vec)

for (l in 2:19) {
  genelist <- panc.gene[which(panc.panc > quant[l] & panc.panc <= quant[l+1])]

  glA <- geneB[which(geneA %in% genelist)] 
  glB <- geneA[which(geneB %in% genelist)]
  gl.all <- c(glA, glB)
  gene.count <- count(gl.all)
  orfs <- gene.count$x
  freqs <- gene.count$freq
  vec <- numeric(length=cc.dim)

#Note that this loop is necessary to avoid miss counting owing to the redundancy of genes
  for (i in 1:length(orfs)) {
    orf.ith <- orfs[i]
    orf.freq <- freqs[i]
    orf.term <- cc.GO.term[which(cc.GO.gene %in% orf.ith)]
    for (j in 1:length(orf.term)) {
      na <- which(cc.go.cat %in% orf.term[j])
      vec[na] <- vec[na] + orf.freq
    } 
  }
all.count <- c(all.count, vec)
}

  genelist.20 <- panc.gene[which(panc.panc > quant[20])]
  glA <- geneB[which(geneA %in% genelist.20)] 
  glB <- geneA[which(geneB %in% genelist.20)]
  gl.all <- c(glA, glB)
  gene.count <- count(gl.all)
  orfs <- gene.count$x
  freqs <- gene.count$freq
  vec <- numeric(length=cc.dim)

#Note that this loop is necessary to avoid miss counting owing to the redundancy of genes
  for (i in 1:length(orfs)) {
    orf.ith <- orfs[i]
    orf.freq <- freqs[i]
    orf.term <- cc.GO.term[which(cc.GO.gene %in% orf.ith)]
    for (j in 1:length(orf.term)) {
      na <- which(cc.go.cat %in% orf.term[j])
      vec[na] <- vec[na] + orf.freq
    } 
  }
  all.count <- c(all.count, vec)

A <- matrix(all.count, nrow = cc.dim, ncol = 20)

write.table(A, file="human.rnf43.all.cc.txt", sep=",", col.names=F, row.names=F, quote=F)

# Now the ms02star permutations
for (p in 1:100) {
permutation.file <- paste("../ms02star/human/", "ms02.", p, ".csv", sep="")
permutation <- read.csv(permutation.file, header=T, stringsAsFactors = F)
geneA <- permutation$id1
geneB <- permutation$id2

perm.count <- c()

#gene lists Q1 to Q20
  glA <- geneB[which(geneA %in% genelist.1)] 
  glB <- geneA[which(geneB %in% genelist.1)]
  gl.all <- c(glA, glB)
  gene.count <- count(gl.all)
  orfs <- gene.count$x
  freqs <- gene.count$freq
  vec <- numeric(length=cc.dim)

#Note that this loop is necessary to avoid miss counting owing to the redundancy of genes
  for (i in 1:length(orfs)) {
    orf.ith <- orfs[i]
    orf.freq <- freqs[i]
    orf.term <- cc.GO.term[which(cc.GO.gene %in% orf.ith)]
    for (j in 1:length(orf.term)) {
      na <- which(cc.go.cat %in% orf.term[j])
      vec[na] <- vec[na] + orf.freq
    } 
  }
  perm.count <- c(perm.count, vec)

for (l in 2:19) {
  genelist <- panc.gene[which(panc.panc > quant[l] & panc.panc <= quant[l+1])]
  glA <- geneB[which(geneA %in% genelist)] 
  glB <- geneA[which(geneB %in% genelist)]
  gl.all <- c(glA, glB)
  gene.count <- count(gl.all)
  orfs <- gene.count$x
  freqs <- gene.count$freq
  vec <- numeric(length=cc.dim)

#Note that this loop is necessary to avoid miss counting owing to the redundancy of genes
  for (i in 1:length(orfs)) {
    orf.ith <- orfs[i]
    orf.freq <- freqs[i]
    orf.term <- cc.GO.term[which(cc.GO.gene %in% orf.ith)]
    for (j in 1:length(orf.term)) {
      na <- which(cc.go.cat %in% orf.term[j])
      vec[na] <- vec[na] + orf.freq
    } 
  }
perm.count <- c(perm.count, vec)
}

  glA <- geneB[which(geneA %in% genelist.20)] 
  glB <- geneA[which(geneB %in% genelist.20)]
  gl.all <- c(glA, glB)
  gene.count <- count(gl.all)
  orfs <- gene.count$x
  freqs <- gene.count$freq
  vec <- numeric(length=cc.dim)

#Note that this loop is necessary to avoid miss counting owing to the redundancy of genes
  for (i in 1:length(orfs)) {
    orf.ith <- orfs[i]
    orf.freq <- freqs[i]
    orf.term <- cc.GO.term[which(cc.GO.gene %in% orf.ith)]
    for (j in 1:length(orf.term)) {
      na <- which(cc.go.cat %in% orf.term[j])
      vec[na] <- vec[na] + orf.freq
    } 
  }
  perm.count <- c(perm.count, vec)

B <- matrix(perm.count, nrow = cc.dim, ncol = 20)

output <- paste("ms02.human", "/", "rnf43.heatmap", "/", "ms02.", p, ".cc.matrix.csv", sep="")

write.table(B, file=output, sep=",", col.names=F, row.names=F, quote=F)
}

### this block perform z-score calculations
library("microbenchmark")
library("matrixStats")

conn.dim <- 20
hspin <- matrix(as.numeric(unlist(read.table("human.rnf43.all.cc.txt", header=F, sep=","))), nrow=cc.dim, ncol=conn.dim)
obs <- c(hspin)

perm <- c()
for (i in 1:100) {
    name <- paste("ms02.human", "/", "rnf43.heatmap", "/", "ms02.", i, ".cc.matrix.csv", sep="")
    mat <- matrix(as.numeric(unlist(read.table(name, header=F, sep=","))), nrow=cc.dim, ncol=conn.dim)
    perm <- rbind(perm, c(mat))
}

mean <- colMeans(perm)
std <- colSds(perm)

zscore <- round((obs - mean)/std, 3)

z <- matrix(zscore, nrow=cc.dim, ncol=conn.dim)

write.table(z, file="human.rnf43.cc.NP.z.csv", sep=",", row.names=F, col.names=F, quote=F)

library('gplots')
z <- t(z)

rownames(z) <- c("Q1","Q2","Q3","Q4","Q5","Q6","Q7","Q8","Q9","Q10",
                 "Q11","Q12","Q13","Q14","Q15","Q16","Q17","Q18","Q19","Q20")

colnames(z) <- cc.go.cat

colors = c(seq(min(z),-10.1,length=100),seq(-9.9,9.9,length=100),seq(10.1,max(z),length=100))
#colors = c(seq(min(z), max(z), length=300))
my_palette <- colorRampPalette(c("blue2", "white", "red2"))(n = 299)

png(filename = "human.rnf43.cc.NP.png",width=6, height=5.5, res=1200, unit="in")
heatmap.2(z, col=my_palette, breaks=colors, Rowv=F,
          trace='none', offsetRow = 0, offsetCol = 0,
          xlab="Biological Process Terms", ylab="Quantiles of Bayesian Factors",
          margins = c(2,3.5), key.title = "Color Bar", key.xlab="Z-score", key.ylab=NA,
          revC = T,
          labCol = NA, #labRow =,
          #srtCol=45, adjCol=c(1,0),
          #lmat=rbind(c(0,3,4), c(2,1,0)), lwid=c(1.5,4,2),
          scale="none", dendrogram = "col", symbreaks=T, symm=F, symkey = F)
dev.off()

hm <- heatmap.2(z, col=my_palette, breaks=colors, Rowv=F,
          trace='none', offsetRow = 0, offsetCol = 0,
          xlab="Biological Process Terms", ylab="Quantiles of Bayesian Factors",
          margins = c(2,3.5), key.title = "Color Bar", key.xlab="Z-score", key.ylab=NA,
          revC = T,
          labCol = NA, 
          scale="none", dendrogram = "col", symbreaks=T, symm=F, symkey = F)
hc.col <- as.hclust(hm$colDendrogram)
pdf("human.rnf43.cc.NP.tree.pdf", width=150, height=4,paper='special')
plot(hc.col, xlab="BP Terms", main="Z-scores, Hierachical Clustering", cex=.8)
dev.off()

