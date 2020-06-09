rm(list=ls())
#library(plyr)
#library(gplots)
library(microbenchmark)
library(matrixStats)

#all GO terms
random.dim <- 50

int.matrix <- matrix(as.numeric(unlist(read.table("hhi.csv", header=F, sep=","))), ncol=1)

int.perm <- c()
for (i in 1:10000) {
    int.name <- paste("ms02.hhi/", "ms02.", i, ".hhi.csv", sep="")
    int.mat <- matrix(as.numeric(unlist(read.table(int.name, header=F, sep=","))), ncol=1)
    int.perm <- rbind(int.perm, c(int.mat))
}

int.mean <- colMeans(int.perm)
int.std <- colSds(int.perm)

int.z <- (int.matrix - int.mean)/int.std

z.matrix <- matrix(int.z, nrow=random.dim, ncol=random.dim)
mean.matrix <- matrix(int.mean, nrow=random.dim, ncol=random.dim)
std.matrix <- matrix(int.std, nrow=random.dim, ncol=random.dim)

write.table(z.matrix, file="hhi.z.csv", sep=",", row.names = F, col.names = F, quote=F)
write.table(mean.matrix, file="hhi.mean.csv", sep=",", row.names = F, col.names = F, quote=F)
write.table(std.matrix, file="hhi.sd.csv", sep=",", row.names = F, col.names = F, quote=F)

