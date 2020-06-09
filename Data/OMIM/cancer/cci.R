#
library(gplots) #for heatmap.2

cancer <- read.csv("cancer.gene.csv", header=T, stringsAsFactors = F)
list <- unique(cancer$NAME[!is.na(cancer$SN)])
cancer.dim <- length(list)

cancer.abb <- c("AE", "BL", "BR", "CO", "CR", "HE",
                "ES", "GA", "TY", "LE", "OV", "PR",
                "SL", "TH", "LU", "PA", "MI", "CE",
                "BO", "EN", "OL", "SE")

cci.dat <- read.csv("cci.z.csv", header=F, stringsAsFactors = F)
cci.z <- matrix(unlist(cci.dat), nrow=cancer.dim, ncol=cancer.dim)

colnames(cci.z) <- list
rownames(cci.z) <- list

my_palette <- colorRampPalette(c("blue2", "white", "red2"))(n = 20)
#colors= c(seq(-5,-1.5, length=10), seq(-1.49, 1.49, length=10), seq(1.5,5, length=10))
colors = c(seq(-10,10,length=21))

png("cci.z.heatmap.png", width=12, height=11, res=600, units="in")
heatmap.2(cci.z, col=my_palette, trace='none', breaks=colors, 
          key.xlab=NA, key.title="Z-score", key.ylab=NA, key.xtickfun = NULL, key.ytickfun = NULL,
          adjCol=c(1,0.5), dendrogram = "both", srtCol=45, 
          margins=c(14,18.5), sepwidth=c(0.01,0.01), symbreaks = TRUE,
          sepcolor="grey", colsep=1:cancer.dim, rowsep=1:cancer.dim)
dev.off()

cci.z.mat <- as.matrix(cci.z)
cci.z.net <- graph.adjacency(cci.z.mat, mode="undirected", weighted=T, diag=F)

top5 <- quantile(cci.z[which(!is.na(cci.z))], probs=seq(0,1,1/10))[10]
bottom5 <- quantile(cci.z[which(!is.na(cci.z))], probs=seq(0,1,1/10))[2]

E(cci.z.net)$color <- ifelse(E(cci.z.net)$weight > 0, "red", "blue")
coloring <- E(cci.z.net)$color
cci.weight <- ifelse((E(cci.z.net)$weight > top5 | (E(cci.z.net)$weight < bottom5)), 
                     abs(E(cci.z.net)$weight), -0.5)

pdf("cci.z.network.pdf", height=8, width=8, paper='special')
plot.igraph(cci.z.net, vertex.label=cancer.abb, layout=layout_in_circle,
            edge.color = coloring, edge.width=cci.weight)
dev.off()
