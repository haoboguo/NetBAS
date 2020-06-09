#
library(gplots) #for heatmap.2
library(igraph)

#all diseases (8753)
list.file <- read.csv("all.ids.csv", header=T, stringsAsFactors = F)
list <- list.file$ID[which(list.file$GN > 14)]

disorder.dim <- length(list)

#print the list
#refined.disorders <- c("ID", "Name", "GN")
#for (i in 1:disorder.dim) {
#list.name <- list.file$Name[which(list.file$ID == list[i])]
#list.gn <- list.file$GN[which(list.file$ID == list[i])]
#refined.disorders <- rbind(refined.disorders,
#                           c(list[i], list.name, list.gn))
#}
#write.table(refined.disorders, file="refined.disorders.csv",
#            col.names = F, row.names = F, sep=",", quote=T)

#diseases with >=15 genes (645)
disorder.summary <- read.csv("refined.disorders.csv", header=T,
                             stringsAsFactors = F)


#calculating overlaps between disorders
disorder.gene <- read.csv("Discurate.full.non.MIR.csv", header=T, stringsAsFactors = F)

#overlap.matrix <- matrix(0, nrow=disorder.dim, ncol=disorder.dim)
#id.list <- disorder.summary$ID
#for (i in 1:length(id.list)) {
#  gene.list.a <- disorder.gene$gene[which(disorder.gene$ID %in% id.list[i])]
#  for (j in 1:length(id.list)) {
#    gene.list.b <- disorder.gene$gene[which(disorder.gene$ID %in% id.list[j])]
#    overlap <- length(which(gene.list.a %in% gene.list.b))
#    overlap.matrix[i,j] = overlap.matrix[i,j] + overlap/(sqrt(length(gene.list.a)*length(gene.list.b)))
#  }
#}

#row.names(overlap.matrix) <- disorder.summary$Name
#colnames(overlap.matrix) <- disorder.summary$Name

#my_palette <- colorRampPalette(c("blue2", "white", "red2"))(n = 49)
#colors = c(seq(0,0.25,length=50))

#png("overlap.png", width=12, height=12, units = "in", res=600)
#heatmap.2(overlap.matrix, col=my_palette, trace='none', breaks=colors, 
#          key.xlab=NA, key.title="Normalized Overlap", key.ylab=NA, key.xtickfun = NULL, key.ytickfun = NULL,
#          srtCol=45, adjCol=c(1,0), dendrogram = "both",
#          margins=c(14,18.5))
#dev.off()

#write.table(overlap.matrix, file="Overlap.matrix.csv", col.names = F, row.names = F, quote=T, sep=",")

#length(which(overlap.matrix == 0)) #229430 among 416025 disorder pairs, 55.14813% 
#length(which(overlap.matrix == 1)) #1861, there are certain redundant disorders

overlap.matrix <- read.csv("Overlap.matrix.csv", header=F, stringsAsFactors = F)
colnames(overlap.matrix) <- disorder.summary$Name
row.names(overlap.matrix) <- disorder.summary$Name

##Jaccard indices
#jac.matrix <- matrix(0, nrow=disorder.dim, ncol=disorder.dim)
#id.list <- disorder.summary$ID
#for (i in 1:length(id.list)) {
#  gene.list.a <- disorder.gene$gene[which(disorder.gene$ID %in% id.list[i])]
#  for (j in 1:length(id.list)) {
#    gene.list.b <- disorder.gene$gene[which(disorder.gene$ID %in% id.list[j])]
#    a.and.b <- length(which(gene.list.a %in% gene.list.b))
#    a.or.b  <- length(unique(c(gene.list.b, gene.list.a)))
#    jac.matrix[i,j] = jac.matrix[i,j] + a.and.b/a.or.b
#  }
#}

#row.names(jac.matrix) <- disorder.summary$Name
#colnames(jac.matrix) <- disorder.summary$Name

#my_palette <- colorRampPalette(c("blue2", "white", "red2"))(n = 49)
#colors = c(seq(0,0.1,length=50))

#png("jaccard.png", width=12, height=12, units = "in", res=600)
#heatmap.2(jac.matrix, col=my_palette, trace='none', breaks=colors, 
#          key.xlab=NA, key.title="Jaccard Indices", key.ylab=NA, key.xtickfun = NULL, key.ytickfun = NULL,
#          srtCol=45, adjCol=c(1,0), dendrogram = "both",
#          margins=c(14,18.5))
#dev.off()

#write.table(jac.matrix, file="Jaccard.indices.csv", col.names = F, row.names = F, quote=T, sep=",")

#length(which(jac.matrix == 0)) #229430 among 416025 disorder pairs, 55.14813% 
#length(which(jac.matrix == 1)) #1861, there are certain redundant disorders
##!the same as overlap matrix above


####disorder genes
disorder.gene.uniq <- unique(disorder.gene$gene[which(disorder.gene$ID %in% disorder.summary$ID)])
# 5894 uniq genes in all 645 disorders

#gene.counts <- c("gene", "disorder.number")
#for (i in 1:length(disorder.gene.uniq)) {
#  gene <- disorder.gene.uniq[i]
#  number <- 0
#  for (j in 1:length(disorder.summary$ID)) {
#    genes.list <- disorder.gene$gene[which(disorder.gene$ID %in% disorder.summary$ID[j])]
#    number <- number + length(which(genes.list %in% gene))
#  }
#  gene.counts <- rbind(gene.counts, c(gene, number))
#}
#write.table(gene.counts, file="gene.counts.csv", sep="\t", col.names = F, row.names = F, quote=T)

gene.counts <- read.csv("gene.counts.csv", sep="\t", header=T, stringsAsFactors = F)
pdf("gene.counts.histogram.pdf", width=5, height=4, paper='special')
hist(gene.counts$disorder.number, breaks=100, xlab="Disorder Number")
dev.off()

gene.counts$gene[which(gene.counts$disorder.number >100)]
#"IL1B"  "IL6"   "SOD2"  "TNF"   "TP53"  "PTGS2"
gene.counts$gene[which(gene.counts$disorder.number >50)]
# "APOE"   "BDNF"   "CNR1"   "ACE"    "DRD2"   "ESR1"   "FGFR1"  "IL1B"   "IL6"    "MMP9"  
# "MTHFR"  "NPY"    "SOD2"   "TNF"    "TP53"   "BCL2"   "HMOX1"  "IGF1"   "INS"    "LEP"   
# "NOS3"   "PPARG"  "VEGFA"  "APC"    "BRAF"   "CTNNB1" "EGFR"   "ERBB2"  "IFNG"   "KRAS"  
# "MYC"    "NOS2"   "PTGS2"  "STAT3"  "AGT"    "NGF"    "POMC"   "SOD1"   "ALB"    "IFNA2" 
# "TGFB1"  "CCL2"   "MET"    "PIK3CA" "PTEN"   "FOS"    "CAT"    "CSF3"   "CSF2"
length(which(gene.counts$disorder.number == 1))
# 1791 genes only belong to one disorder; likewise, 1075 genes belong to two disorders
# and 650 genes to three disorders, ...

###the z-scores
ddi.dat <- read.csv("ddi.z.1k.csv", header=F, stringsAsFactors = F)
ddi.z <- matrix(unlist(ddi.dat), nrow=disorder.dim, ncol=disorder.dim)

colnames(ddi.z) <- disorder.summary$Name
row.names(ddi.z) <- disorder.summary$Name

my_palette <- colorRampPalette(c("blue2", "white", "red2"))(n = 20)
#colors= c(seq(-5,-1.5, length=10), seq(-1.49, 1.49, length=10), seq(1.5,5, length=10))
colors = c(seq(-10,10,length=21))

hm <- heatmap.2(ddi.z, col=my_palette, trace='none', breaks=colors, 
                key.xlab=NA, key.title="Z-score", key.ylab=NA, key.xtickfun = NULL, key.ytickfun = NULL,
                srtCol=45, adjCol=c(1,0), dendrogram = "both",
                margins=c(14,18.5))
hc <- as.hclust(hm$rowDendrogram)
pdf("disorder.tree.pdf", height=10, width=80, paper='special')
plot(hc, xlab="Disorder ID", cex=.8)
dev.off()

hierarchy.order <- hc$order
disorder.summary$Name[hierarchy.order[[1]]]

####disorder network
ddi.z.mat <- as.matrix(ddi.z)
ddi.z.net <- graph.adjacency(ddi.z.mat, mode="undirected", weighted=T, diag=F)

top1 <- quantile(ddi.z, probs=seq(0,1,1/100))[100]
bottom1 <- quantile(ddi.z, probs=seq(0,1,1/50))[2]

E(ddi.z.net)$color <- ifelse(E(ddi.z.net)$weight > 0, "red", "blue")
coloring <- E(ddi.z.net)$color
E(ddi.z.net)$weight <- ifelse((E(ddi.z.net)$weight > 20), 
                     abs(E(ddi.z.net)$weight), -0.5)
#ddi.weight <- abs(E(ddi.z.net)$weight)

pdf("ddi.z.network.pdf", height=8, width=8, paper='special')
plot.igraph(ddi.z.net, vertex.label=NA, layout=layout_in_circle,
            edge.color = coloring, edge.width=E(ddi.z.net)$weight/10, vertex.size=4)
dev.off()

ddi.z.plus <- ddi.z

ddi.z.mat.plus <- as.matrix(ddi.z.plus)
ddi.z.net.plus <- graph.adjacency(ddi.z.mat.plus, mode="undirected", weighted=T, diag=F)
gc <- cluster_fast_greedy(ddi.z.net.plus, weights=E(ddi.z.net.plus)$weight, merges=T,
                          modularity = T, membership = T)

wc <- cluster_walktrap(ddi.z.net.plus, weights=E(ddi.z.net.plus)$weight, merges=T,
                       modularity = T, membership=T)

sc <- cluster_spinglass(ddi.z.net.plus, weights=E(ddi.z.net.plus)$weight, vertex=NULL,
                        spins=25, implementation = c("neg"), gamma=0.1, gamma.minus=0.1)

disgenet.table <- read.csv("../all_gene_disease_associations.tsv", sep="\t", header=T,
                           stringsAsFactors = F)
#disease.types <- unique(disgenet.table$diseaseType)
#disease.class <- unique(disgenet.table$diseaseClass)
disease.semantictype    <- unique(disgenet.table$diseaseSemanticType)
#the disease semantic type (25 types) will be used to classify the diseases

#name.type <- c("Disorder.Name", "Semantic.Type", "Disease.class")
#for (i in 1:disorder.dim) {
#  disorder.name <- disorder.summary$Name[i]
#  semantictype <- disgenet.table$diseaseSemanticType[which(disgenet.table$diseaseName == disorder.name)[1]]
#  disclass <- disgenet.table$diseaseClass[which(disgenet.table$diseaseName == disorder.name)[1]]
#  name.type <- rbind(name.type, c(disorder.name, semantictype, disclass))
#}
#write.table(name.type, file="Dis.Name.v.Type-class.csv", col.names = F, row.names = F, quote=T, sep=",")

library(Corbi)

###Cancers
name.type <- read.csv("Dis.Name.v.Type-class.csv", header=T, stringsAsFactors = F)
cancers <- which(name.type$Semantic.Type == "Neoplastic Process")
cancer.z.matrix <- submatrix(ddi.z, rows=cancers, cols=cancers)
colnames(cancer.z.matrix) <- disorder.summary$Name[cancers]
row.names(cancer.z.matrix) <- disorder.summary$Name[cancers]

my_palette <- colorRampPalette(c("blue2", "white", "red2"))(n = 20)
#colors= c(seq(-5,-1.5, length=10), seq(-1.49, 1.49, length=10), seq(1.5,5, length=10))
colors = c(seq(-10,10,length=21))
png("cancers.png", width=13, height=12, units = "in", res=600)
heatmap.2(cancer.z.matrix, col=my_palette, trace='none', breaks=colors, 
          key.xlab=NA, key.title="Z-scores", key.ylab=NA, key.xtickfun = NULL, key.ytickfun = NULL,#          srtCol=45, adjCol=c(1,0), dendrogram = "both",
          margins=c(14,18.5), srtCol=45, adjCol=c(1,0), dendrogram = "both")
         # sepcolor="grey", colsep=1:length(cancers), rowsep=1:length(cancers))
dev.off()

hm2 <- heatmap.2(cancer.z.matrix, col=my_palette, trace='none', breaks=colors, 
                key.xlab=NA, key.title="Z-score", key.ylab=NA, key.xtickfun = NULL, key.ytickfun = NULL,
                srtCol=45, adjCol=c(1,0), dendrogram = "both",
                margins=c(14,18.5))
hc2 <- as.hclust(hm2$rowDendrogram)

pdf("cancer.tree.pdf", height=10, width=30, paper='special')
plot(hc2, xlab="Disorders", cex=.8)
dev.off()

cancer.branch1 <- hc2$order[1:42]
cancer.branch1.matrix <- submatrix(cancer.z.matrix, rows=cancer.branch1, cols = cancer.branch1)
colnames(cancer.branch1.matrix) <- hc2$labels[cancer.branch1]
row.names(cancer.branch1.matrix) <- hc2$labels[cancer.branch1]
png("cancer.branch1.png", width=13, height=12, units = "in", res=600)
heatmap.2(cancer.branch1.matrix, col=my_palette, trace='none', breaks=colors, 
          key.xlab=NA, key.title="Z-scores", key.ylab=NA, key.xtickfun = NULL, key.ytickfun = NULL,#          srtCol=45, adjCol=c(1,0), dendrogram = "both",
          margins=c(14,18.5), srtCol=45, adjCol=c(1,0), dendrogram = "both",
          sepcolor="grey", colsep=1:42, rowsep=1:42)
dev.off()

##Mental disorders

mental <- which(name.type$Semantic.Type == "Mental or Behavioral Dysfunction")
mental.z.matrix <- submatrix(ddi.z, rows=mental, cols=mental)
colnames(mental.z.matrix) <- disorder.summary$Name[mental]
row.names(mental.z.matrix) <- disorder.summary$Name[mental]

png("mental.png", width=13, height=12, units = "in", res=600)
heatmap.2(mental.z.matrix, col=my_palette, trace='none', breaks=colors, 
          key.xlab=NA, key.title="Z-scores", key.ylab=NA, key.xtickfun = NULL, key.ytickfun = NULL,#          srtCol=45, adjCol=c(1,0), dendrogram = "both",
          margins=c(14,18.5), srtCol=45, adjCol=c(1,0), dendrogram = "both")
dev.off()

hm3 <- heatmap.2(mental.z.matrix, col=my_palette, trace='none', breaks=colors, 
                 key.xlab=NA, key.title="Z-score", key.ylab=NA, key.xtickfun = NULL, key.ytickfun = NULL,
                 srtCol=45, adjCol=c(1,0), dendrogram = "both",
                 margins=c(14,18.5))
hc3 <- as.hclust(hm3$rowDendrogram)

pdf("mental.tree.pdf", height=10, width=30, paper='special')
plot(hc3, xlab="Disorders", cex=.8)
dev.off()

mental.branch2 <- hc3$order[22:68]
mental.branch2.matrix <- submatrix(mental.z.matrix, rows=mental.branch2, cols = mental.branch2)
colnames(mental.branch2.matrix) <- hc3$labels[mental.branch2]
row.names(mental.branch2.matrix) <- hc3$labels[mental.branch2]
png("mental.branch2.png", width=13, height=12, units = "in", res=600)
heatmap.2(mental.branch2.matrix, col=my_palette, trace='none', breaks=colors, 
          key.xlab=NA, key.title="Z-scores", key.ylab=NA, key.xtickfun = NULL, key.ytickfun = NULL,#          srtCol=45, adjCol=c(1,0), dendrogram = "both",
          margins=c(14,18.5), srtCol=45, adjCol=c(1,0), dendrogram = "both",
          sepcolor="grey", colsep=1:47, rowsep=1:47)
dev.off()


##Congenital Abnormality
Congen <- which(name.type$Semantic.Type == "Congenital Abnormality")
Congen.z.matrix <- submatrix(ddi.z, rows=Congen, cols=Congen)
colnames(Congen.z.matrix) <- disorder.summary$Name[Congen]
row.names(Congen.z.matrix) <- disorder.summary$Name[Congen]
png("Congenital.png", width=12, height=12, units = "in", res=600)
heatmap.2(Congen.z.matrix, col=my_palette, trace='none', breaks=colors, 
          key.xlab=NA, key.title="Z-scores", key.ylab=NA, key.xtickfun = NULL, key.ytickfun = NULL,#          srtCol=45, adjCol=c(1,0), dendrogram = "both",
          margins=c(14,18.5), srtCol=45, adjCol=c(1,0), dendrogram = "both")
dev.off()

##Other
other <- which(name.type$Semantic.Type == "Disease or Syndrome")
other.z.matrix <- submatrix(ddi.z, rows=other, cols=other)
colnames(other.z.matrix) <- disorder.summary$Name[other]
row.names(other.z.matrix) <- disorder.summary$Name[other]
png("other.png", width=12, height=12, units = "in", res=600)
heatmap.2(other.z.matrix, col=my_palette, trace='none', breaks=colors, 
          key.xlab=NA, key.title="Z-scores", key.ylab=NA, key.xtickfun = NULL, key.ytickfun = NULL,#          srtCol=45, adjCol=c(1,0), dendrogram = "both",
          margins=c(14,18.5), srtCol=45, adjCol=c(1,0), dendrogram = "both")
dev.off()


##C04-C06, cancer in digestive system
c04.list <- which(name.type$Disease.class=="C04;C06" |
                  name.type$Disease.class=="C04;C06;C16" |
                  name.type$Disease.class=="C04;C06;C16;C18" |
                  name.type$Disease.class=="C04;C06;C19")
c04.z.matrix <- submatrix(ddi.z, rows=c04.list, cols=c04.list)
colnames(c04.z.matrix) <- disorder.summary$Name[c04.list]
row.names(c04.z.matrix) <- disorder.summary$Name[c04.list]
png("c04.c06.png", width=13, height=12, units = "in", res=600)
heatmap.2(c04.z.matrix, col=my_palette, trace='none', breaks=colors, 
          key.xlab=NA, key.title="Z-scores", key.ylab=NA, key.xtickfun = NULL, key.ytickfun = NULL,#          srtCol=45, adjCol=c(1,0), dendrogram = "both",
          margins=c(14,18.5), srtCol=45, adjCol=c(1,0), dendrogram = "both",
          sepcolor="grey", colsep=1:42, rowsep=1:42)
dev.off()

name.type$Disorder.Name[which(name.type$Disease.class=="C04;C15;C20")]

##C4-C15-C20-C10
c15.list <- which(name.type$Disease.class=="C04;C15" |
                  name.type$Disease.class=="C04;C15;C20" |
                  name.type$Disease.class=="C02;C04;C15;C20" |
                  name.type$Disease.class=="C04;C14;C15;C20" |
                  name.type$Disease.class=="C04;C10" |
                  name.type$Disease.class=="C04;C10;C16")

c15.z.matrix <- submatrix(ddi.z, rows=c15.list, cols=c15.list)
colnames(c15.z.matrix) <- disorder.summary$Name[c15.list]
row.names(c15.z.matrix) <- disorder.summary$Name[c15.list]
png("c15.c20.c10.png", width=13.5, height=12, units = "in", res=600)
heatmap.2(c15.z.matrix, col=my_palette, trace='none', breaks=colors, 
          key.xlab=NA, key.title="Z-scores", key.ylab=NA, key.xtickfun = NULL, key.ytickfun = NULL,#          srtCol=45, adjCol=c(1,0), dendrogram = "both",
          margins=c(16,22), srtCol=45, adjCol=c(1,0), dendrogram = "both",
          sepcolor="grey", colsep=1:42, rowsep=1:42)
dev.off()

##C4
c4.only <- which(name.type$Disease.class=="C04")
c4.only.z.matrix <- submatrix(ddi.z, rows=c4.only, cols=c4.only)
colnames(c4.only.z.matrix) <- disorder.summary$Name[c4.only]
row.names(c4.only.z.matrix) <- disorder.summary$Name[c4.only]

########
## we'll measure two specific sets:
##     C0009402 for "Colorectal Carcinoma"
## vs C0678222 for "Breast Carcinoma"

pin <- read.csv("../../human.pin.csv", header=T, stringsAsFactors = F)

geneA <- pin$geneA
geneB <- pin$geneB

id.a <- c("C0009402")  #colorectal
id.b <- c("C0678222")  #breast
gene.list.a <- disorder.gene$gene[which(disorder.gene$ID %in% id.a)]
gene.list.b <- disorder.gene$gene[which(disorder.gene$ID %in% id.b)]

gene.list.a[which(gene.list.a %in% gene.list.b)]

#"AKT1", "CASP8", "CDH1", "CTNNB1", "EP300", "KRAS", "MMP1", "PIK3CA", "TP53", "EXO1", "CHEK2"
# Is it able to identify the driver gene or distinguish the ONG vs TRG?
ab.int <- which(((geneA %in% gene.list.a) & (geneB %in% gene.list.b)) |
                     ((geneA %in% gene.list.b) & (geneB %in% gene.list.a)))
subA <- geneA[ab.int]
subB <- geneB[ab.int]

subnet <- data.frame(cbind(subA, subB))
sub.graph <- graph.data.frame(subnet, directed=F)

'%ni%' <- Negate('%in%')

blue.id <- which((as_ids(V(sub.graph)) %in% gene.list.b) &
                   (as_ids(V(sub.graph)) %ni% gene.list.a))
red.id <-  which((as_ids(V(sub.graph)) %in% gene.list.a) &
                   (as_ids(V(sub.graph)) %ni% gene.list.b))
yellow.id <- which((as_ids(V(sub.graph)) %in% gene.list.b) &
                     (as_ids(V(sub.graph)) %in% gene.list.a))
sub.order <- V(sub.graph)[c(blue.id, yellow.id, red.id)]
coords <- layout_in_circle(sub.graph, order = sub.order)

color <- rep("NA", times=length(V(sub.graph)))
color[red.id] <- rep("red", times=length(red.id))
color[blue.id] <- rep("blue", times=length(blue.id))
color[yellow.id] <- rep("yellow", times=length(yellow.id))
V(sub.graph)$color <- color

pdf("Breast-Colorectal.subnetwork.pdf")
plot.igraph(sub.graph, vertex.color=V(sub.graph)$color,
            vertex.size=8, edge.width=1,vertex.label=NA,
            order=sub.order,layout=coords)
dev.off()
length(E(sub.graph))  #=407


######## ms02
ms02 <- read.csv("../../../ms02star/human/ms02.1.csv", header=T, stringsAsFactors = F)

ms02.geneA <- ms02$id1
ms02.geneB <- ms02$id2

ms02.ab.int <- which(((ms02.geneA %in% gene.list.a) & (ms02.geneB %in% gene.list.b)) |
                  ((ms02.geneA %in% gene.list.b) & (ms02.geneB %in% gene.list.a)))
ms02.subA <- ms02.geneA[ms02.ab.int]
ms02.subB <- ms02.geneB[ms02.ab.int]

ms02.subnet <- data.frame(cbind(ms02.subA, ms02.subB))
ms02.sub.graph <- graph.data.frame(ms02.subnet, directed=F)

'%ni%' <- Negate('%in%')

ms02.blue.id <- which((as_ids(V(ms02.sub.graph)) %in% gene.list.b) &
                   (as_ids(V(ms02.sub.graph)) %ni% gene.list.a))
ms02.red.id <-  which((as_ids(V(ms02.sub.graph)) %in% gene.list.a) &
                   (as_ids(V(ms02.sub.graph)) %ni% gene.list.b))
ms02.yellow.id <- which((as_ids(V(ms02.sub.graph)) %in% gene.list.b) &
                     (as_ids(V(ms02.sub.graph)) %in% gene.list.a))
ms02.sub.order <- V(ms02.sub.graph)[c(ms02.blue.id, ms02.yellow.id, ms02.red.id)]
ms02.coords <- layout_in_circle(ms02.sub.graph, order = ms02.sub.order)

ms02.color <- rep("NA", times=length(V(ms02.sub.graph)))
ms02.color[ms02.red.id] <- rep("red", times=length(ms02.red.id))
ms02.color[ms02.blue.id] <- rep("blue", times=length(ms02.blue.id))
ms02.color[ms02.yellow.id] <- rep("yellow", times=length(ms02.yellow.id))
V(ms02.sub.graph)$color <- ms02.color

pdf("ms02.Breast-Colorectal.subnetwork.pdf")
plot.igraph(ms02.sub.graph, vertex.color=V(ms02.sub.graph)$color,
            vertex.size=8, edge.width=1,vertex.label=NA,
            order=sub.order,layout=coords)
dev.off()
length(E(ms02.sub.graph))  #=213

########
#barplot test
c4.c6.nodiag <- c4.z.matrix
diag(c4.c6.nodiag) <- NA
pdf("C4.C6.z.boxplot.pdf", width=10, height=6, paper='special')
par(mar=c(12,3,1,1))
boxplot(c4.c6.nodiag, main="Neoplasms (C04) + Digestive System Diseases (C06)", xaxt="n", col='brown')
text(seq(1,22), par("usr")[3]-0.25, srt=45, adj=1, xpd=T, 
     col='blue',labels=paste(rownames(c4.c6.nodiag)), cex=0.8, line=10)
dev.off()




