# Venn diagram for Biogrid/CellMap (negative and positive interactions are treated separatedly

library(VennDiagram)

david.file <- read.csv("david.bp.id.txt", header=F)
david.id <- david.file$V1[1:200]
net.file <- read.csv("HMK.bp.full.txt", sep="\t", header=T)
net.id <- net.file$GO.ID[1:200]
web.file <- read.csv("WebGestalt.2.txt", sep="\t", header=T)
web.id <- web.file$goId[1:200]

area1 <- length(net.id)
area2 <- length(david.id)
area3 <- length(web.id)
net.david <- net.id[which(net.id %in% david.id)]
net.web <- net.id[which(net.id %in% web.id)]
david.web <- david.id[which(david.id %in% web.id)]
n12 <- length(net.david)
n13 <- length(net.web)
n23 <- length(david.web)
net.david.web <- net.id[which(net.id %in% david.web)]
n123 <- length(net.david.web)

venn.plot.triple <- draw.triple.venn(area1, area2, area3, n12, n23, n13, n123,
                                     category = c("NetBAS", "DAVID", "WebGestAlt"),
                                     col=c("red","blue","green"), fill=c("red","blue","green"),
                                     euler.d=T, scaled=T,
                                     cat.col=c("red","blue","green"),
                                     cat.fontface=2,
                                     sep.dist=0.1, rotation.degree=30)

pdf("bp200.venn.3way.pdf", width=3.5, height=3.5, paper='special')
grid.draw(venn.plot.triple)
dev.off()

