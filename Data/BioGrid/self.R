my.data <- read.csv("yeast.biogrid.csv", header=T)
geneA <- my.data$geneA
geneB <- my.data$geneB
pair<-c("geneA","geneB")
for (i in 1:length(geneA)) {
  orfa <- as.character(geneA[i])
  orfb <- as.character(geneB[i])
  if (orfa == orfb) {
     pA <- NULL
     pB <- NULL
  } else {
     pA <- orfa
     pB <- orfb
  }
  pair <- rbind(pair,sort(c(pA,pB)))
}
write.table(pair,file="yeast.self.csv", sep=",", row.names=F, col.names=F, quote=T)
