hset.file <- "h.all.v6.2.symbols.gmt.txt"
hset.Line <- file(hset.file, open="r")
hset.data <- readLines(hset.Line)
for (i in 1:length(hset.data)) {
  single.set <- unlist(strsplit(hset.data[i], "\t"))
  file.name <- paste(single.set[1], ".csv", sep="")
  genes <- c("gene")
  for (j in 3:length(single.set)) {
    genes <- rbind(genes, single.set[j])
  }
  write.table(genes, file=file.name, col.names = F, row.names = F, quote=T)
}