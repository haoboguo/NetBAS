##

bp <- read.csv("DAVID.bp.all.txt", sep="\t", header=T, stringsAsFactors=F)
cc <- read.csv("DAVID.cc.all.txt", sep="\t", header=T, stringsAsFactors=F)
mf <- read.csv("DAVID.mf.all.txt", sep="\t", header=T, stringsAsFactors=F)

bp.new <- cbind(bp$Term[which(bp$PValue < 0.05)], bp$PValue[which(bp$PValue < 0.05)])
cc.new <- cbind(cc$Term[which(cc$PValue < 0.05)], cc$PValue[which(cc$PValue < 0.05)])
mf.new <- cbind(mf$Term[which(mf$PValue < 0.05)], mf$PValue[which(mf$PValue < 0.05)])

write.table(bp.new, file="david.bp.new.txt", sep="\t", col.names=F, row.names=F, quote=F)
write.table(cc.new, file="david.cc.new.txt", sep="\t", col.names=F, row.names=F, quote=F)
write.table(mf.new, file="david.mf.new.txt", sep="\t", col.names=F, row.names=F, quote=F)
