# Venn diagram for Biogrid/CellMap (negative and positive interactions are treated separatedly

library(VennDiagram)

david.file <- read.csv("../david.bp.id.txt", header=F)
david.id <- david.file$V1
net.file <- read.csv("../HMK.bp.full.txt", sep="\t", header=T)
net.id <- net.file$GO.ID[which(net.file$Z.score > 10)]
web.file <- read.csv("HMK.bp.full.txt", sep="\t", header=T)
web.id <- web.file$GO.ID[which(web.file$Z.score > 10)]
perm.file <- read.csv("david.bp.id.txt", header=F)
perm.id <- perm.file$V1

area1 <- length(net.id)
area2 <- length(david.id)
area3 <- length(web.id)
area4 <- length(perm.id)
net.david <- net.id[which(net.id %in% david.id)]
net.web <- net.id[which(net.id %in% web.id)]
net.perm <- net.id[which(net.id %in% perm.id)]
david.web <- david.id[which(david.id %in% web.id)]
david.perm <- david.id[which(david.id %in% perm.id)]
web.perm <- web.id[which(web.id %in% perm.id)]
n12 <- length(net.david)
n13 <- length(net.web)
n14 <- length(net.perm)
n23 <- length(david.web)
n24 <- length(david.perm)
n34 <- length(web.perm)
net.david.web <- net.id[which(net.id %in% david.web)]
n123 <- length(net.david.web)
n124 <- length(net.id[which(net.id %in% david.perm)])
n134 <- length(net.id[which(net.id %in% web.perm)])
n234 <- length(david.id[which(david.id %in% web.perm)])
n1234 <- length(perm.id[which(perm.id %in% net.david.web)])

venn.plot.quadraple <- draw.quad.venn(area1, area2, area3, area4,
                                      n12, n13, n14, n23, n24, n34,
                                      n123, n124, n134, n234, n1234,
                                      category = c("NetBAS", "DAVID", "netbas.c3","david.c3"),
                                      fill = c("red", "blue", "green", "orange"),
#                                      cex = 1, cat.cex = 1,
                                      cat.fontface=2,
                                      cat.col = c("red", "blue", "green", "orange"))

pdf("bp.venn.4way.2.pdf", width=3.5, height=3.5, paper='special')
grid.draw(venn.plot.quadraple)
dev.off()
