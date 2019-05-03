# Venn diagram for Biogrid/CellMap (negative and positive interactions are treated separatedly

library(VennDiagram)

david.file <- read.csv("david.path.csv", header=T)
david.id <- david.file$ID
net.file <- read.csv("netbas.path.csv", header=T)
net.id <- net.file$ID
path.file <- read.csv("pathfindR.path.csv", header=T)
path.id <- path.file$ID

area1 <- length(net.id)
area2 <- length(david.id)
area3 <- length(path.id)
net.david <- net.id[which(net.id %in% david.id)]
net.path <- net.id[which(net.id %in% path.id)]
david.path <- david.id[which(david.id %in% path.id)]
n12 <- length(net.david)
n13 <- length(net.path)
n23 <- length(david.path)
net.david.path <- net.id[which(net.id %in% david.path)]
n123 <- length(net.david.path)

venn.plot.triple <- draw.triple.venn(area1, area2, area3, n12, n23, n13, n123,
                                     category = c("NetBAS", "DAVID", "pathfindR"),
                                     col=c("red","blue","green"), fill=c("red","blue","green"),
                                     euler.d=T, scaled=T,
                                     cat.col=c("red","blue","green"),
                                     cat.fontface=2,
                                     sep.dist=0.1, rotation.degree=30)

pdf("kegg.venn.3way.pdf", width=3.5, height=3.5, paper='special')
grid.draw(venn.plot.triple)
dev.off()

