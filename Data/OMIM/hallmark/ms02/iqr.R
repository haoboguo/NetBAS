rm(list=ls())
#library(plyr)
#library(gplots)
library(microbenchmark)
library(matrixStats)

#all GO terms
random.dim <- 50

int.matrix <- as.numeric(unlist(read.table("hhi.csv", header=F, sep=",")))

int.perm <- c()
for (i in 1:10000) {
    int.name <- paste("ms02.hhi/", "ms02.", i, ".hhi.csv", sep="")
    int.mat <- matrix(as.numeric(unlist(read.table(int.name, header=F, sep=","))), ncol=1)
    int.perm <- rbind(int.perm, c(int.mat))
}

iqr.median <- c()
iqr.del <- c()
iqr.score <- c()
for (j in 1:2500) {
  int.dat <- int.perm[,j]
  int.quantile <- quantile(int.dat)
  int.delta <- int.quantile[4] - int.quantile[2]
  int.iqr <- (int.matrix[j] - median(int.dat)) / int.delta
  iqr.median <- rbind(iqr.median, median(int.dat))
  iqr.del <- rbind(iqr.del, int.delta)
  iqr.score <- rbind(iqr.score, int.iqr)
}

iqr.matrix <- matrix(iqr.score, nrow=random.dim, ncol=random.dim)
median.matrix <- matrix(iqr.median, nrow=random.dim, ncol=random.dim)
delta.matrix <- matrix(iqr.del, nrow=random.dim, ncol=random.dim)

write.table(iqr.matrix, file="hhi.iqr.csv", sep=",", row.names = F, col.names = F, quote=F)
write.table(median.matrix, file="hhi.median.csv", sep=",", row.names = F, col.names = F, quote=F)
write.table(delta.matrix, file="hhi.delta.csv", sep=",", row.names = F, col.names = F, quote=F)

