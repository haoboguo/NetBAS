# Venn diagram for Biogrid/CellMap (negative and positive interactions are treated separatedly

library(VennDiagram)

david.file <- read.csv("david.bp.id.txt", header=F)
david.id <- david.file$V1[1:100]
net.file <- read.csv("HMK.bp.new.txt", sep="\t", header=T)
net.id <- net.file$GO.ID[1:100]
web.file <- read.csv("WebGestalt.2.txt", sep="\t", header=T)
web.id <- web.file$goId[1:100]
perm.file <- read.csv("HMK.perm.bp.full.txt", sep="\t", header=T)
perm.id <- perm.file$GO.ID[1:100]

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

venn.plot <- draw.pairwise.venn(area1,area2,cross.area=n12,
             category=c("NetBAS", "DAVID"),
             col=c("red","blue"), fill=c("red","blue"),
             cat.pos = c(0, 180),
             euler.d = TRUE, sep.dist = 0.03,rotation.degree = 0)

pdf("bp.venn.2way.pdf", width=3, height=3, paper='special') 
grid.draw(venn.plot)
dev.off()

venn.plot.triple <- draw.triple.venn(area1, area2, area3, n12, n23, n13, n123,
                                     category = c("NetBAS", "DAVID", "WebGestAlt"),
                                     col=c("red","blue","yellow"), fill=c("red","blue","yellow"),
                                     euler.d=T, scaled=T,
                                     #cat.col=c("red","blue","yellow"),
                                     cat.fontface=2,
                                     sep.dist=0.1, rotation.degree=30)

pdf("bp.venn.3way.pdf", width=3, height=3, paper='special')
grid.draw(venn.plot.triple)
dev.off()

venn.plot.quadraple <- draw.quad.venn(area1, area2, area3, area4,
                                      n12, n13, n14, n23, n24, n34,
                                      n123, n124, n134, n234, n1234,
                                      category = c("NetBAS", "DAVID", "WebGestAlt","NodeRep"),
                                      fill = c("red", "blue", "yellow", "grey"),
                                      cex = 1, cat.cex = 1,
                                      cat.col = c("red", "blue", "yellow", "grey"))

pdf("bp.venn.4way.pdf", width=3.5, height=3.5, paper='special')
grid.draw(venn.plot.quadraple)
dev.off()
