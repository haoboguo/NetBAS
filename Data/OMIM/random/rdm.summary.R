rri1.100 <- matrix(unlist(read.csv("random1-100/rri.z.csv", header=F, stringsAsFactors = F)), nrow=50)
rri100 <- matrix(unlist(read.csv("random100/rri.z.csv", header=F, stringsAsFactors = F)), nrow=50)
rri15 <- matrix(unlist(read.csv("random15/rri.z.csv", header=F, stringsAsFactors = F)), nrow=50)
rri50 <- matrix(unlist(read.csv("random50/rri.z.csv", header=F, stringsAsFactors = F)), nrow=50)
rri200 <- matrix(unlist(read.csv("random200/rri.z.csv", header=F, stringsAsFactors = F)), nrow=50)
rri15.200 <- matrix(unlist(read.csv("random15-200/rri.z.csv", header=F, stringsAsFactors = F)), nrow=50)

diag(rri1.100) = NA
diag(rri100) = NA
diag(rri15) = NA
diag(rri200) = NA
diag(rri15.200) = NA

rri.rdm <- c(rri15, rri50, rri200, rri100, rri15.200, rri1.100)

length(rri.rdm)

length(which(rri.rdm > 2)) / length(rri.rdm) # 0.030

length(which(rri.rdm < -2)) / length(rri.rdm) # 0.011

pdf("rdm.z.hist.pdf", height=4.5, width=6, paper='special')
hist(rri.rdm, breaks=50, xlab="Z-score", main="", xlim=c(-4,4), col=c("green"))
dev.off()

rri.rdm.nodiag <- rri.rdm[which(!is.na(rri.rdm))]
length(rri.rdm.nodiag) #14750
length(which(rri.rdm.nodiag > 2))/length(rri.rdm.nodiag) #0.030
length(which(rri.rdm.nodiag < -2))/length(rri.rdm.nodiag) #0.011
