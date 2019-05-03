#Table1

tab <- read.csv("table1b.csv", header=T, stringsAsFactors=F)
#Set,V,E,Iso,Dmax,CC,Clust,Cliq,R,B,Y,JBP,JCC,JMF,BPmax,CCmax,MFmax
ccoe <- (tabel1$CC-min(tabel1$CC))/(max(tabel1$CC) - min(tabel1$CC))
clust <- (tabel1$Clust-min(tabel1$Clust))/(max(tabel1$Clust) - min(tabel1$Clust))
clq <- (tabel1$Cliq-min(tabel1$Cliq))/(max(tabel1$Cliq) - min(tabel1$Cliq))
bp <- (tabel1$BPmax-min(tabel1$BPmax))/(max(tabel1$BPmax) - min(tabel1$BPmax))
cc <- (tabel1$CCmax-min(tabel1$CCmax))/(max(tabel1$CCmax) - min(tabel1$CCmax))
mf <- (tabel1$MFmax-min(tabel1$MFmax))/(max(tabel1$MFmax) - min(tabel1$MFmax))
jbp <- (tabel1$JBP-min(tabel1$JBP))/(max(tabel1$JBP) - min(tabel1$JBP))
jcc <- (tabel1$JCC-min(tabel1$JCC))/(max(tabel1$JCC) - min(tabel1$JCC))
jmf <- (tabel1$JMF-min(tabel1$JMF))/(max(tabel1$JMF) - min(tabel1$JMF))

info <- cbind(ccoe, clust, clq, bp, cc, mf, jbp, jcc, jmf)
colnames(info) <- c("Clust.Coeff", "Max.Clust", "Max.Clique", "BPmax", "CCmax", "MFmax", "Jbp", "Jcc", "Jmf")

pairs(info)
