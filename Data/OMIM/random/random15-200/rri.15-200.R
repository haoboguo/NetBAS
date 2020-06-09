#
library(gplots) #for heatmap.2

list <- c(1:50)
random.dim <- 50

rri.dat <- read.csv("rri.z.csv", header=F, stringsAsFactors = F)
rri.z <- matrix(unlist(rri.dat), nrow=random.dim, ncol=random.dim)

colnames(rri.z) <- list
rownames(rri.z) <- list

my_palette <- colorRampPalette(c("blue2", "white", "red2"))(n = 20)
#colors= c(seq(-5,-1.5, length=10), seq(-1.49, 1.49, length=10), seq(1.5,5, length=10))
colors = c(seq(-10,10,length=21))

png("rri15-200.z.heatmap.new3.png", width=12, height=11, res=600, units="in")
heatmap.2(rri.z, col=my_palette, trace='none', breaks=colors, 
          key.xlab=NA, key.title="Interaction Z-score Heatmap", key.ylab=NA, key.xtickfun = NULL, key.ytickfun = NULL,
          #srtCol=45, adjCol=c(1,0), dendrogram = "both",
          Rowv = F, Colv = F, revC=T,
          margins=c(14,18.5), sepwidth=c(0.01,0.01), symbreaks = TRUE,
          sepcolor="grey", colsep=1:random.dim, rowsep=1:random.dim)
dev.off()

rri.z.mat <- as.matrix(rri.z)
diag(rri.z.mat) <- NA
rri.z.net <- graph.adjacency(rri.z.mat, mode="undirected", weighted=T, diag=F)

top2 <- quantile(rri.z, probs=seq(0,1,1/50))[50]
bottom2 <- quantile(rri.z, probs=seq(0,1,1/50))[2]

E(rri.z.net)$color <- ifelse(E(rri.z.net)$weight > 0, "red", "blue")
coloring <- E(rri.z.net)$color
rri.weight <- ifelse((E(rri.z.net)$weight > 2 | (E(rri.z.net)$weight < -2)), 
                     abs(E(rri.z.net)$weight), -0.5)

size.file <- read.csv("set.size.csv", header=T)
sizes <- size.file$size

pdf("rri15-200.z.network.pdf", height=8, width=8, paper='special')
plot.igraph(rri.z.net, vertex.label=c(1:50), layout=layout_in_circle,
            edge.color = coloring, edge.width=rri.weight/4, vertex.size=sqrt(sizes-1)*1.2)
dev.off()

rri.nodiag <- rri.z
diag(rri.nodiag) <- NA
boxplot(rri.nodiag)
pdf("rri.z.boxplot.nodiag.pdf", height=6, width=8, paper='special')
par(mar=c(15,5,1,1))
boxplot(rri.nodiag, main="Random-Random Z-score BoxPlot", xaxt="n", col='brown')
text(seq(1,50), par("usr")[3]-0.5, srt=-45, adj=1, xpd=T, 
     labels=c(1:50), cex=0.75)
dev.off()