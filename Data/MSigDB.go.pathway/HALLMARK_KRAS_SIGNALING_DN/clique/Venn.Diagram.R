# Venn diagram for Biogrid/CellMap (negative and positive interactions are treated separatedly

library(VennDiagram)

david.file <- read.csv("david.bp.id.txt", header=F)
david.id <- david.file$V1
web.file <- read.csv("WebGestalt.txt", sep="\t", header=T)
web.id <- web.file$goId
net.file <- read.csv("HMK.bp.new.txt", sep="\t", header=T)
net.id <- net.file$GO.ID[1:100]

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

venn.plot <- draw.pairwise.venn(area1=area1,area2=area2,cross.area=n12,
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
                                     sep.dist=0.1)

pdf("bp.venn.3way.pdf", width=3.5, height=3.5, paper='special')
grid.draw(venn.plot.triple)
dev.off()

