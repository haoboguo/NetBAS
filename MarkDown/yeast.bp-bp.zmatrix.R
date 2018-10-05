library("microbenchmark")
library("matrixStats")
GOcategory.file <- read.csv("../Data/yeast.bp.term.csv",header=TRUE, stringsAsFactors=F)
go.cat <- GOcategory.file$GO.term
dim <- length(go.cat)

biogrid <- matrix(as.numeric(unlist(read.table("yeast.bp.matrix.csv", header = F, sep = ","))), ncol=dim)
obs <- c(biogrid)

perm <- c()
for (i in 1:100) {
    name <- paste("ms02.",  i, "bp", ".matrix.csv", sep = "")
    mat <- matrix(as.numeric(unlist(read.table(name, header=F, sep=","))), ncol=dim)
    perm <- rbind(perm, c(mat))
}

mean <- colMeans(perm)
std <- colSds(perm)

zscore <- round((obs - mean)/std, 3)

z <- t(matrix(zscore, ncol=dim))

write.table(z, file="yeast.bp-bp.zscore.100.csv", sep = ",", row.names=F, col.names=F)


