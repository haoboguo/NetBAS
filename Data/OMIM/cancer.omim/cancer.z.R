#
library(igraph)
library(gplots) #for heatmap.2
cancer <- read.csv("cancer.gene.csv", header=T, stringsAsFactors = F)

cancer.name <- unique(cancer$NAME[!is.na(cancer$SN)])

#calculate association between a pair of two cancers

cancer.dim <- length(cancer.name)

cci.dat <- read.csv("cci.z.csv", header=F, stringsAsFactors = F)
cci.z <- matrix(unlist(cci.dat), nrow=cancer.dim, ncol=cancer.dim)


colnames(cci.z) <- cancer.name
rownames(cci.z) <- cancer.name

my_palette <- colorRampPalette(c("blue2", "white", "red2"))(n = 20)
colors= c(seq(-8,8, length=21))

png("cci.z.heatmap.new2.png", width=12, height=4, res=600, units="in")
heatmap.2(cci.z, col=my_palette, trace='none', breaks=colors, 
          key.xlab=NA, key.title="Interaction Z-score Heatmap", key.ylab=NA, key.xtickfun = NULL, key.ytickfun = NULL,
          srtCol=45, adjCol=c(1,0), dendrogram = "both",
          margins=c(14,18.5), sepwidth=c(0.01,0.01), symbreaks = TRUE,
          sepcolor="grey", colsep=1:length(cancer.name), rowsep=1:length(cancer.name))
dev.off()

cci.z.mat <- as.matrix(cci.z)
cci.z.net <- graph.adjacency(cci.z.mat, mode="undirected", weighted=T, diag=F)
E(cci.z.net)$color <- ifelse(E(cci.z.net)$weight > 0, "brown", "green")

coloring <- E(cci.z.net)$color
cci.weight <- ifelse(abs(E(cci.z.net)$weight) > 2, abs(E(cci.z.net)$weight), -0.5)

# in two-tailed test Z > 1.96 corresponds to p < 0.05
# Z > 1.6449 is for p < 0.10
# Z > 1.2819 for p < 0.10 in one-tailed test
# Z > 1.6449 for p < 0.50 in one-tailed test

radi <- c()
for (i in 1:length(cancer.name)) {
  radi <- rbind(radi, sqrt(length(which(cancer$NAME %in% cancer.name[i]))))
}

black.id <- c(1:13)
blue.id <- c(14:16)
red.id <- c(18,17,19)

sub.order <- V(cci.z.net)[c(red.id, black.id, blue.id)]
coords <- layout_in_circle(cci.z.net, order = sub.order)

color <- rep("NA", times=19)
color[red.id] <- rep("red", times=length(red.id))
color[blue.id] <- rep("orange", times=length(blue.id))
color[black.id] <- rep("grey", times=length(black.id))
V(cci.z.net)$color <- color

pdf("cci.net.pdf", width=8, height=8, paper='special')
plot.igraph(cci.z.net, layout=coords, cex=0.6,
            edge.color = coloring, edge.width=cci.weight*2, vertex.size=radi*6)
dev.off()

cci.z.nodiag <- cci.z
diag(cci.z.nodiag) <- NA

pdf("cci.z.boxplot.pdf", height=4, width=6, paper='special')
par(mar=c(8,3,1,1))
boxplot(cci.z.nodiag, main="Z-score Box Plot", 
        space=0, xaxt="n", col='brown')
text(seq(1,19), par("usr")[3]-0.25, srt=45, adj=1, xpd=T, 
     col='blue',labels=cancer.name, cex=1)
dev.off()


