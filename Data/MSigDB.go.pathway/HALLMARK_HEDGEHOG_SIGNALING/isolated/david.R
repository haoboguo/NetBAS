##

bp <- read.csv("DAVID.bp.all.txt", sep="\t", header=T, stringsAsFactors=F)

bp.new <- cbind(bp$Term[which(bp$PValue < 0.05)], bp$PValue[which(bp$PValue < 0.05)])

write.table(bp.new, file="david.bp.new.txt", sep="\t", col.names=F, row.names=F, quote=F)
