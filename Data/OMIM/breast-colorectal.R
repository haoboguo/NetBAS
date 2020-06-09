#library(gplots) #for heatmap.2
library(igraph)
cancer <- read.csv("cancer.gene.csv", header=T, stringsAsFactors = F)

cancer.name <- unique(cancer$NAME[!is.na(cancer$SN)])

pin <- read.csv("../human.pin.csv", header=T, stringsAsFactors = F)

geneA <- pin$geneA
geneB <- pin$geneB

breast.list <- cancer$Gene[which(cancer$NAME == cancer.name[3])]
colorectal.list <- cancer$Gene[which(cancer$NAME == cancer.name[5])]

breast.list[which(breast.list %in% colorectal.list)]
#"AKT1"   "PIK3CA" "TP53"

br.co.int <- which(((geneA %in% breast.list) & (geneB %in% colorectal.list)) | 
                  ((geneA %in% colorectal.list) & (geneB %in% breast.list)))

subA <- geneA[br.co.int]
subB <- geneB[br.co.int]

subnet <- data.frame(cbind(subA, subB))
sub.graph <- graph.data.frame(subnet, directed=F)
#plot.igraph(sub.graph)

'%ni%' <- Negate('%in%')

blue.id <- which((as_ids(V(sub.graph)) %in% breast.list) & 
                 (as_ids(V(sub.graph)) %ni% colorectal.list))
red.id <-  which((as_ids(V(sub.graph)) %in% colorectal.list) & 
                  (as_ids(V(sub.graph)) %ni% breast.list))
yellow.id <- which((as_ids(V(sub.graph)) %in% colorectal.list) & 
                     (as_ids(V(sub.graph)) %in% breast.list))
sub.order <- V(sub.graph)[c(blue.id, yellow.id, red.id)]
coords <- layout_in_circle(sub.graph, order = sub.order)

color <- rep("NA", times=length(V(sub.graph)))
color[red.id] <- rep("red", times=length(red.id))
color[blue.id] <- rep("blue", times=length(blue.id))
color[yellow.id] <- rep("yellow", times=length(yellow.id))
V(sub.graph)$color <- color


pdf("Breast-Colorectal.subnetwork.pdf")
plot.igraph(sub.graph, vertex.color=V(sub.graph)$color, 
            vertex.size=20, edge.width=2,
            order=sub.order, layout=coords)
dev.off()


library(VennDiagram)
pdf("Breast-Colorecta-Venn.pdf", height=3, width=3, paper='special')
venn.plot <- draw.pairwise.venn(32,21,3,
                                category=c("Colorectal Cancer", "Breast Cancer"),
                                col=c("red","blue"), fill=c("red","blue"),
                                cat.pos = c(0, 0),cat.col=c("red", "blue"),
                                euler.d = TRUE, sep.dist = 0.03,rotation.degree = 0)
dev.off()