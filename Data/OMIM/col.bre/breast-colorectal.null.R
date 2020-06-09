#library(gplots) #for heatmap.2
library(igraph)
cancer <- read.csv("cancer.gene.csv", header=T, stringsAsFactors = F)

cancer.name <- unique(cancer$NAME[!is.na(cancer$SN)])

pin <- read.csv("../../../ms02star/human/ms02.1.csv", header=T, stringsAsFactors = F)

geneA <- pin$id1
geneB <- pin$id2

breast.list <- cancer$Gene[which(cancer$NAME == cancer.name[1])]
colorectal.list <- cancer$Gene[which(cancer$NAME == cancer.name[2])]

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

breast.sub <- as_ids(V(sub.graph))[c(blue.id, yellow.id)]
breast.iso <- breast.list[which(breast.list %ni% breast.sub)]
colorectal.sub <- as_ids(V(sub.graph))[c(red.id, yellow.id)]
colorectal.iso <- colorectal.list[which(colorectal.list %ni% colorectal.sub)]
blue.id <- c(blue.id, 25:29)
red.id <- c(red.id, 30:45)

sub.graph = add_vertices(sub.graph, length(breast.iso), name=breast.iso)
sub.graph = add_vertices(sub.graph, length(colorectal.iso), name=colorectal.iso)

sub.order <- V(sub.graph)[c(blue.id, yellow.id, red.id)]
coords <- layout_in_circle(sub.graph, order = sub.order)

color <- rep("NA", times=length(V(sub.graph)))
color[red.id] <- rep("red", times=length(red.id))
color[blue.id] <- rep("blue", times=length(blue.id))
color[yellow.id] <- rep("yellow", times=length(yellow.id))
V(sub.graph)$color <- color


pdf("Breast-Colorectal.null.subnetwork.pdf")
plot.igraph(sub.graph, vertex.color=V(sub.graph)$color, 
            vertex.size=15, edge.width=1, cex=0.6,
            order=sub.order, layout=coords)
dev.off()


