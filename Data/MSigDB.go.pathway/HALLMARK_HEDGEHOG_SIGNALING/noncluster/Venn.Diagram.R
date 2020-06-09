# Venn diagram for Biogrid/CellMap (negative and positive interactions are treated separatedly

library(VennDiagram)

david.file <- read.csv("david.bp.new.txt", sep="\t", header=T, stringsAsFactors=F)
david.id <- david.file$ID[which(david.file$pValue < 0.01)]
web.file <- read.csv("WebGestalt.txt", sep="\t", header=T, stringsAsFactors=F)
web.id <- web.file$goId[which(web.file$pValue < 0.01)]
net.file <- read.csv("HMK.bp.full.txt", sep="\t", header=T, stringsAsFactors=F)
net.id <- net.file$GO.ID[which(net.file$Z.score > 5)]

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
             category=c("NetPAS", "DAVID"),
             col=c("red","blue"), fill=c("red","blue"),
             cat.pos = c(0, 180),
             euler.d = TRUE, sep.dist = 0.03,rotation.degree = 0)

pdf("bp.venn.2way.pdf", width=3, height=3, paper='special') 
grid.draw(venn.plot)
dev.off()

venn.plot.triple <- draw.triple.venn(area1, area2, area3, n12, n23, n13, n123,
                                     category = c("NetPAS", "DAVID", "WebGestAlt"),
                                     col=c("red","blue","yellow"), fill=c("red","blue","yellow"),
                                     euler.d=T, scaled=T,
                                    # cat.col=c("red","blue","yellow"),
                                     cat.fontface=2,
                                     sep.dist=0.1)

pdf("bp.venn.3way.pdf", width=3, height=3, paper='special')
grid.draw(venn.plot.triple)
dev.off()

