#
library(gplots) #for heatmap.2

list.file <- read.csv("../../../MSigDB.go.pathway/list", header=F, stringsAsFactors = F)
list <- list.file$V1
hallmark.dim <- 50

hhi.dat <- read.csv("hhi.z.csv", header=F, stringsAsFactors = F)
hhi.z <- matrix(unlist(hhi.dat), nrow=hallmark.dim, ncol=hallmark.dim)

my_palette <- colorRampPalette(c("blue2", "white", "red2"))(n = 20)
colors = c(seq(-10,10,length=21))

hallmark.matrix <- hhi.z
colnames(hallmark.matrix) <- c(1:50)
row.names(hallmark.matrix) <- c(1:50)
png("hhi.z.png", width=12, height=11, res=600, units="in")
heatmap.2(hallmark.matrix, col=my_palette, trace='none', breaks=colors, 
          Rowv = F, Colv = F, revC = T,
          key.xlab=NA, key.title="Association Z-score", 
          key.ylab=NA, key.xtickfun = NULL, key.ytickfun = NULL,
          #srtCol=45, adjCol=c(1,0), 
          dendrogram = "row",
          margins=c(14,18.5), sepwidth=c(0.01,0.01), #symbreaks = TRUE,
          sepcolor="grey", colsep=1:hallmark.dim, rowsep=1:hallmark.dim)
#keysize=1, key.par=list(mar=c(2,1,2,2)))
#lmat=rbind( c(0, 3, 4), c(2,1,1.5)), lwid=c(3, 4, 2))
#lmat=rbind(c(5, 4, 2), c(6, 1, 3)), lhei=c(2, 5), lwid=c(1, 10, 1))
dev.off()

hhi.iqr.dat <- read.csv("hhi.iqr.csv", header=F, stringsAsFactors = F)
hhi.iqr <- matrix(unlist(hhi.iqr.dat), nrow=hallmark.dim, ncol=hallmark.dim)

my_palette <- colorRampPalette(c("blue2", "white", "red2"))(n = 20)
colors = c(seq(-10,10,length=21))

iqr.matrix <- hhi.iqr
colnames(iqr.matrix) <- c(1:50)
row.names(iqr.matrix) <- c(1:50)
png("hhi.iqr.png", width=12, height=11, res=600, units="in")
heatmap.2(iqr.matrix, col=my_palette, trace='none', breaks=colors, 
          Rowv = F, Colv = F, revC = T,
          key.xlab=NA, key.title="Modified Z-score", 
          key.ylab=NA, key.xtickfun = NULL, key.ytickfun = NULL,
          #srtCol=45, adjCol=c(1,0), 
          dendrogram = "row",
          margins=c(14,18.5), sepwidth=c(0.01,0.01), #symbreaks = TRUE,
          sepcolor="grey", colsep=1:hallmark.dim, rowsep=1:hallmark.dim)
#keysize=1, key.par=list(mar=c(2,1,2,2)))
#lmat=rbind( c(0, 3, 4), c(2,1,1.5)), lwid=c(3, 4, 2))
#lmat=rbind(c(5, 4, 2), c(6, 1, 3)), lhei=c(2, 5), lwid=c(1, 10, 1))
dev.off()

png("hhi.iqr.colorbar.png", width=12, height=4, res=600, units="in")
heatmap.2(iqr.matrix, col=my_palette, trace='none', breaks=colors, 
          Rowv = F, Colv = F, revC = T,
          key.xlab=NA, key.title="Modified Z-score", 
          key.ylab=NA, key.xtickfun = NULL, key.ytickfun = NULL,
          #srtCol=45, adjCol=c(1,0), 
          dendrogram = "row",
          margins=c(14,18.5), sepwidth=c(0.01,0.01), #symbreaks = TRUE,
          sepcolor="grey", colsep=1:hallmark.dim, rowsep=1:hallmark.dim)
#keysize=1, key.par=list(mar=c(2,1,2,2)))
#lmat=rbind( c(0, 3, 4), c(2,1,1.5)), lwid=c(3, 4, 2))
#lmat=rbind(c(5, 4, 2), c(6, 1, 3)), lhei=c(2, 5), lwid=c(1, 10, 1))
dev.off()

