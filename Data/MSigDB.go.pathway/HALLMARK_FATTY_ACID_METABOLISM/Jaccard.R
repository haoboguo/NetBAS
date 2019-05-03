# Venn diagram for Biogrid/CellMap (negative and positive interactions are treated separatedly

library(VennDiagram)

jaccard.index <- c("BP.J", "CC.J", "MF.J")

david.bp <- read.csv("david.bp.id.txt", header=F)
david.bp.id <- david.bp$V1
net.bp <- read.csv("HMK.bp.new.txt", sep="\t", header=T)
net.bp.id <- net.bp$GO.ID[which(net.bp$zscore > 5.0)]
david.cc <- read.csv("david.cc.id.txt", header=F)
david.cc.id <- david.cc$V1
net.cc <- read.csv("HMK.cc.new.txt", sep="\t", header=T)
net.cc.id <- net.cc$GO.ID[which(net.cc$zscore > 5.0)]
david.mf <- read.csv("david.mf.id.txt", header=F)
david.mf.id <- david.mf$V1
net.mf <- read.csv("HMK.mf.new.txt", sep="\t", header=T)
net.mf.id <- net.mf$GO.ID[which(net.mf$zscore > 5.0)]


bpj <- length(which(net.bp.id %in% david.bp.id))/length(unique(sort(c(net.bp.id, david.bp.id))))
ccj <- length(which(net.cc.id %in% david.cc.id))/length(unique(sort(c(net.cc.id, david.cc.id))))
mfj <- length(which(net.mf.id %in% david.mf.id))/length(unique(sort(c(net.mf.id, david.mf.id))))

jaccard.index <- rbind(jaccard.index,c(bpj, ccj, mfj))

print(jaccard.index)
