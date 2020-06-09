library(MASS)

## fitting a normal distribution
fitnormal <- function (x, exact = TRUE) {
  if (exact) {
    ################################################
    ## Exact inference based on likelihood theory ##
    ################################################
    ## minimum negative log-likelihood (maximum log-likelihood) estimator of `mu` and `phi = sigma ^ 2`
    n <- length(x)
    mu <- sum(x) / n
    phi <- crossprod(x - mu)[1L] / n  # (a bised estimator, though)
    ## inverse of Fisher information matrix evaluated at MLE
    invI <- matrix(c(phi, 0, 0, phi * phi), 2L,
                   dimnames = list(c("mu", "sigma2"), c("mu", "sigma2")))
    ## log-likelihood at MLE
    loglik <- -(n / 2) * (log(2 * pi * phi) + 1)
    ## return
    return(list(theta = c(mu = mu, sigma2 = phi), vcov = invI, loglik = loglik, n = n))
    }
  else {
    ##################################################################
    ## Numerical optimization by minimizing negative log-likelihood ##
    ##################################################################
    ## negative log-likelihood function
    ## define `theta = c(mu, phi)` in order to use `optim`
    nllik <- function (theta, x) {
      (length(x) / 2) * log(2 * pi * theta[2]) + crossprod(x - theta[1])[1] / (2 * theta[2])
      }
    ## gradient function (remember to flip the sign when using partial derivative result of log-likelihood)
    ## define `theta = c(mu, phi)` in order to use `optim`
    gradient <- function (theta, x) {
      pl2pmu <- -sum(x - theta[1]) / theta[2]
      pl2pphi <- -crossprod(x - theta[1])[1] / (2 * theta[2] ^ 2) + length(x) / (2 * theta[2])
      c(pl2pmu, pl2pphi)
      }
    ## ask `optim` to return Hessian matrix by `hessian = TRUE`
    ## use "..." part to pass `x` as additional / further argument to "fn" and "gn"
    ## note, we want `phi` as positive so box constraint is used, with "L-BFGS-B" method chosen
    init <- c(sample(x, 1), sample(abs(x) + 0.1, 1))  ## arbitrary valid starting values
    z <- optim(par = init, fn = nllik, gr = gradient, x = x, lower = c(-Inf, 0), method = "L-BFGS-B", hessian = TRUE)
    ## post processing ##
    theta <- z$par
    loglik <- -z$value  ## flip the sign to get log-likelihood
    n <- length(x)
    ## Fisher information matrix (don't flip the sign as this is the Hessian for negative log-likelihood)
    I <- z$hessian / n  ## remember to take average to get mean
    invI <- solve(I, diag(2L))  ## numerical inverse
    dimnames(invI) <- list(c("mu", "sigma2"), c("mu", "sigma2"))
    ## return
    return(list(theta = theta, vcov = invI, loglik = loglik, n = n))
    }
  }

####
int.perm <- c()
for (i in 1:2000) {
    int.name <- paste("ms02.hhi/", "ms02.", i, ".hhi.csv", sep="")
    int.mat <- matrix(as.numeric(unlist(read.table(int.name, header=F, sep=","))), ncol=1)
    int.perm <- rbind(int.perm, c(int.mat))
}

hhi.fit <- fitnormal(int.perm[,1])


pdf("hhi.hist.pdf", width=4, height=4, paper="special")
hist(int.perm[,1], main="Histogram", prob=T, xlab="Interaction Number", br=20)
curve(dnorm(x, hhi.fit$theta[1], sqrt(hhi.fit$theta[2])), add=TRUE, col=2)
dev.off()

pdf("hhi.qqplot.pdf", width=4, height=4, paper="special")
qqnorm(int.perm[,1])
qqline(int.perm[,1], col=2)
dev.off()

hhi.p <- as.numeric(unlist(read.csv(("hhi.p.csv"), header=F, stringsAsFactors = F)))
hhi.q <- qvalue(p = hhi.p)
hhi.p[which(hhi.p == 0)] <- 3.5e-5
q.value <- hhi.q$qvalues
q.value[which(q.value == 0)] <- 3.5e-5

logp <- -log10(hhi.p)
logq <- -log10(q.value)

q.mat <- matrix(logq, nrow=50, ncol=50)
colnames(q.mat) <- c(1:50)
row.names(q.mat) <- c(1:50)
p.mat <- matrix(logp, nrow=50, ncol=50)
colnames(p.mat) <- c(1:50)
row.names(p.mat) <- c(1:50)

my_palette <- colorRampPalette(c("white", "red2"))(n = 20)
colors = c(seq(0,5,length=21))

png("hhi.q.png", width=12, height=11, res=600, units="in")
heatmap.2(q.mat, col=my_palette, trace='none', breaks=colors, 
          Rowv = F, Colv = F, revC = T,
          key.xlab=NA, key.title="log10(q-value)", 
          key.ylab=NA, key.xtickfun = NULL, key.ytickfun = NULL,
          #srtCol=45, adjCol=c(1,0), 
          dendrogram = "row",
          margins=c(14,18.5), sepwidth=c(0.01,0.01), #symbreaks = TRUE,
          sepcolor="grey", colsep=1:hallmark.dim, rowsep=1:hallmark.dim)
#keysize=1, key.par=list(mar=c(2,1,2,2)))
#lmat=rbind( c(0, 3, 4), c(2,1,1.5)), lwid=c(3, 4, 2))
#lmat=rbind(c(5, 4, 2), c(6, 1, 3)), lhei=c(2, 5), lwid=c(1, 10, 1))
dev.off()

png("hhi.p.png", width=12, height=11, res=600, units="in")
heatmap.2(p.mat, col=my_palette, trace='none', breaks=colors, 
          Rowv = F, Colv = F, revC = T,
          key.xlab=NA, key.title="log10(p-value)", 
          key.ylab=NA, key.xtickfun = NULL, key.ytickfun = NULL,
          #srtCol=45, adjCol=c(1,0), 
          dendrogram = "row",
          margins=c(14,18.5), sepwidth=c(0.01,0.01), #symbreaks = TRUE,
          sepcolor="grey", colsep=1:hallmark.dim, rowsep=1:hallmark.dim)
#keysize=1, key.par=list(mar=c(2,1,2,2)))
#lmat=rbind( c(0, 3, 4), c(2,1,1.5)), lwid=c(3, 4, 2))
#lmat=rbind(c(5, 4, 2), c(6, 1, 3)), lhei=c(2, 5), lwid=c(1, 10, 1))
dev.off()

png("hhi.q.colorbar.png", width=12, height=4, res=600, units="in")
heatmap.2(q.mat, col=my_palette, trace='none', breaks=colors, 
          Rowv = F, Colv = F, revC = T,
          key.xlab=NA, key.title="-log10(q-value)", 
          key.ylab=NA, key.xtickfun = NULL, key.ytickfun = NULL,
          #srtCol=45, adjCol=c(1,0), 
          dendrogram = "row",
          margins=c(14,18.5), sepwidth=c(0.01,0.01), #symbreaks = TRUE,
          sepcolor="grey", colsep=1:hallmark.dim, rowsep=1:hallmark.dim)
#keysize=1, key.par=list(mar=c(2,1,2,2)))
#lmat=rbind( c(0, 3, 4), c(2,1,1.5)), lwid=c(3, 4, 2))
#lmat=rbind(c(5, 4, 2), c(6, 1, 3)), lhei=c(2, 5), lwid=c(1, 10, 1))
dev.off()

png("hhi.p.colorbar.png", width=12, height=4, res=600, units="in")
heatmap.2(p.mat, col=my_palette, trace='none', breaks=colors, 
          Rowv = F, Colv = F, revC = T,
          key.xlab=NA, key.title="-log10(p-value)", 
          key.ylab=NA, key.xtickfun = NULL, key.ytickfun = NULL,
          #srtCol=45, adjCol=c(1,0), 
          dendrogram = "row",
          margins=c(14,18.5), sepwidth=c(0.01,0.01), #symbreaks = TRUE,
          sepcolor="grey", colsep=1:hallmark.dim, rowsep=1:hallmark.dim)
#keysize=1, key.par=list(mar=c(2,1,2,2)))
#lmat=rbind( c(0, 3, 4), c(2,1,1.5)), lwid=c(3, 4, 2))
#lmat=rbind(c(5, 4, 2), c(6, 1, 3)), lhei=c(2, 5), lwid=c(1, 10, 1))
dev.off()

