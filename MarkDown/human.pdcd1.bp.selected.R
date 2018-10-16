##

library('gplots')
pd1.bp <- read.csv("human.pdcd1.bp.selected.csv", sep="\t", header=T)
z <- pd1.bp$Z.score[order(pd1.bp$Z.score)]
go <- pd1.bp$GO.ID[order(pd1.bp$Z.score)]
term <- pd1.bp$GO.Term[order(pd1.bp$Z.score)]
z <- matrix(z, ncol=1)
row.names(z) <- term
colnames(z) <- NA

colors = c(seq(0, 40, length=20))
my_palette <- colorRampPalette(c("white", "red2"))(n = 19)

png(filename = "human.pdcd1.bp.selected.png",width=5, height=5, res=1200, unit="in")
heatmap.2(cbind(z,z), trace="none", Colv=NA, Rowv=NA, dendrogram="none",
          col=my_palette, breaks=colors, revC=T, key.title=NA, key.xlab="Z-score", key.ylab=NA,
          labCol="", labRow="", cellnote=cbind(as.character(go),as.character(go)),notecol=1,
          colsep=1:ncol(z), rowsep=1:nrow(z), sepcolor = "lightgrey")
dev.off()

print(term)
