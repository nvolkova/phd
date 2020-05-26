############################################################
## Cosine similarity simulations for signature comparison ##
############################################################

source('~/Desktop/Git/phd/useful_functions.R')
cosineM <- function(X) {
  return(sapply(1:ncol(X), function(i) sapply(1:ncol(X), function(j) cosine(X[,i], X[,j]))))
}

maxs <- sample(1:5000, 1000, replace=T)
X <- sapply(1:1000, function(x) runif(n = 96, min = 0, max=maxs[x]))
X <- apply(X,2,function(y) y/sum(y))
XM <- cosineM(X)
hist(XM[upper.tri(XM)],breaks=100, main='Distribution of similarities for random signatures', xlab='Similarity score')
abline(v = quantile(as.vector(XM[upper.tri(XM)]), 0.95), col='red', lty=2) # (0.80)

library(bayesm)
Dr3 <- sapply(1:1000, function(x) rdirichlet(alpha = rep(3,96)))
Dr3M <- cosineM(Dr3)
hist(Dr3M[upper.tri(Dr3M)],breaks=100, main='Distribution of similarities for random signatures', xlab='Similarity score')
abline(v = quantile(as.vector(Dr3M[upper.tri(Dr3M)]), 0.95), col='red', lty=2) 

# alphas = 1 -> uniform distribution on simplex
# 2D case
Dr2d <- sapply(1:1000, function(x) rdirichlet(alpha = rep(1,2)))
Dr2dM <- cosineM(Dr2d)
hist(Dr2dM[upper.tri(Dr2dM)],breaks=100, main='Distribution of similarities for random signatures, 2D', xlab='Similarity score')
# 96D case - already roughly normal, but centre not exactly in 0.5 (0.51)
Dr1 <- sapply(1:1000, function(x) rdirichlet(alpha = rep(1,96)))
Dr1M <- cosineM(Dr1)
hist(Dr1M[upper.tri(Dr1M)],breaks=100, main='Distribution of similarities for random signatures', xlab='Similarity score')
abline(v = quantile(as.vector(Dr1M[upper.tri(Dr1M)]), 0.95), col='red', lty=2) # (0.61)
# 1000D case - normal, centre very close to 0.5
Dr1000 <- sapply(1:1000, function(x) rdirichlet(alpha = rep(1,1000)))
Dr1000M <- cosineM(Dr1000)
hist(Dr1000M[upper.tri(Dr1000M)],breaks=100, main='Distribution of similarities for random signatures', xlab='Similarity score')
abline(v = quantile(as.vector(Dr1000M[upper.tri(Dr1000M)]), 0.95), col='red', lty=2)


Dr0.5 <- sapply(1:1000, function(x) rdirichlet(alpha = rep(0.5,96)))
Dr0.5M <- cosineM(Dr0.5)
hist(Dr0.5M[upper.tri(Dr0.5M)],breaks=100, main='Distribution of similarities for random signatures', xlab='Similarity score')
abline(v = quantile(as.vector(Dr0.5M[upper.tri(Dr0.5M)]), 0.95), col='red', lty=2)

Dr0.2 <- sapply(1:1000, function(x) rdirichlet(alpha = rep(0.2,96)))
Dr0.2M <- cosineM(Dr0.2)
hist(Dr0.2M[upper.tri(Dr0.2M)],breaks=100, main='Distribution of similarities for random signatures', xlab='Similarity score')
abline(v = quantile(as.vector(Dr0.2M[upper.tri(Dr0.2M)]), 0.95), col='red', lty=2)


pdf('Simulations.pdf', 15,10)
par(mfrow=c(3,2))
maxs <- sample(1:5000, 1000, replace=T)
X <- sapply(1:1000, function(x) runif(n = 10, min = 0, max=maxs[x]))
X <- apply(X,2,function(y) y/sum(y))
XM <- cosineM(X)
hist(XM[upper.tri(XM)],breaks=100, main='Distribution of similarities for random signatures, 10D', xlab='Similarity score')
abline(v = quantile(as.vector(XM[upper.tri(XM)]), 0.95), col='red', lty=2) # (0.80)
maxs <- sample(1:5000, 1000, replace=T)
X <- sapply(1:1000, function(x) runif(n = 96, min = 0, max=maxs[x]))
X <- apply(X,2,function(y) y/sum(y))
XM <- cosineM(X)
hist(XM[upper.tri(XM)],breaks=100, main='Drawing from uniform distribution and L1-normalising, 96D', xlab='Similarity score')
abline(v = quantile(as.vector(XM[upper.tri(XM)]), 0.95), col='red', lty=2) # (0.80)
maxs <- sample(1:5000, 1000, replace=T)
X <- sapply(1:1000, function(x) runif(n = 300, min = 0, max=maxs[x]))
X <- apply(X,2,function(y) y/sum(y))
XM <- cosineM(X)
hist(XM[upper.tri(XM)],breaks=100, main='Distribution of similarities for random signatures, 300D', xlab='Similarity score')
abline(v = quantile(as.vector(XM[upper.tri(XM)]), 0.95), col='red', lty=2) # (0.80)

hist(Dr2dM[upper.tri(Dr2dM)],breaks=100, main='Drawing from Dir(1), 2D', xlab='Similarity score')
hist(Dr1M[upper.tri(Dr1M)],breaks=100, main='Drawing from Dir(1), 96D', xlab='Similarity score')
hist(Dr1000M[upper.tri(Dr1000M)],breaks=100, main='Drawing from Dir(1), 1000D', xlab='Similarity score')
dev.off()

### Distribution of angles in COSMIC signatures

cosmic <- cosineM(cancer_signatures)
hist(cosmic[upper.tri(cosmic)], breaks=20)
summary(cosmic[upper.tri(cosmic)])
pdf('Similarities_in_COSMIC_signatures.pdf', 12,10)
hist(cosmic[upper.tri(cosmic)], prob=T, breaks=20, main='Similarities between COSMIC signatures', xlab = 'Cosine similarity')
lines(density(cosmic[upper.tri(cosmic)], adjust=2), lty="dotted") 

image.plot(cosmic)
for (x in 1:ncol(cosmic))
  for (y in 1:ncol(cosmic))
    text((x-1)/(ncol(cosmic)-1), (y-1)/(ncol(cosmic)-1), sprintf("%0.2f", cosmic[x,y]))

dev.off()


# Now draw a signature from bootstrapped COADSTAD signatures and compare to randomly drawn within CI worm signature
big.sig <- NULL
for (i in 1:100) {
  load(paste0("/Volumes/yoda2/sigs",i*100,".RData"))
  big.sig <- c(big.sig,sigs)
}

mm <- rbind(COAD.mutation.counts,STAD.mutation.counts)
mut_mat = t(mm) + 0.0001
res <- nmf(mut_mat,rank=8,seed=123456,method='brunet')
defsigs <- NMF::basis(res)
defsigs <- defsigs[,c(8,1,2,5,4,7,3,6)]
colnames(defsigs) <- c("Clock-1", "Clock-2", "POLE", "17-like", "MMR-1", "MMR-2", "MMR-3", "MMR-4")

set = list()
for (k in 1:8) {
  set[[k]] <- sapply(1:10000, function(i) 
    big.sig[[i]][,which.max(sapply(1:8, function(j) cosine(big.sig[[i]][,j],defsigs[,k])))] / sum(big.sig[[i]][,which.max(sapply(1:8, function(j) cosine(big.sig[[i]][,j],defsigs[,k])))]))
}

similarity_hist <- function(defsigs,set,k) {
  hist(sapply(1:ncol(set[[k]]), function(s) cosine(defsigs[,k],set[[k]][,s])), breaks=100)
}
similarity_hist(defsigs,set,5)


hist(sapply(1:ncol(set[[5]]), function(s) cosine(learned.sigs.exome[2,-c(35,39,43,47,97:104)],set[[5]][-c(35,39,43,47,97:104),s])), breaks=100)
hist(sapply(1:ncol(set[[5]]), function(s) cosine(learned.sigs.exome[2,-c(97:104)],set[[5]][-c(97:104),s])), breaks=100)
#hist(sapply(1:ncol(set[[5]]), function(s) cosine(learned.sigs.exome[3,-c(35,39,43,47,97:104)],set[[5]][-c(35,39,43,47,97:104),s])), breaks=100)

cosine(learned.sigs.exome[3,-c(35,39,43,47,97:104)],defsigs[-c(35,39,43,47,97:104),5]) # 0.92
cosine(learned.sigs.exome[2,-c(35,39,43,47,97:104)],defsigs[-c(35,39,43,47,97:104),5])

# Now take the worm sigs randomly from their CIs

load('~/Learned_signatures_indel.RData')
learned.sigs <- nmSolve(worm.donor.mut.mat,small.X,maxIter=10000, tol = 1e-06, div.err = 1e-10)
Y = worm.donor.mut.mat
mut_matrix_lower <- matrix(NA,nrow=96,ncol=nrow(learned.sigs),dimnames=list(colnames(learned.sigs)[1:96],row.names(learned.sigs)))
mut_matrix_upper <- matrix(NA,nrow=96,ncol=nrow(learned.sigs),dimnames=list(colnames(learned.sigs)[1:96],row.names(learned.sigs)))
library(MASS)
for (i in 1:96) {
  cov.mat <- poisI(as.matrix(small.X),learned.sigs[,i],Y[,i])
  to.keep <- colnames(cov.mat)[which(diag(cov.mat)!=0)]
  if (length(to.keep)!=0) {
    cov.mat <- cov.mat[to.keep,to.keep]
    SVD <- svd(cov.mat)$d
    if (length(which(SVD<1e-6))>0)
      cov.mat <- ginv(cov.mat)
    if (length(which(SVD<1e-6))==0)
      cov.mat <- solve(cov.mat)
    row.names(cov.mat) <- to.keep
    colnames(cov.mat) <- to.keep
    sterr <- sqrt(diag(cov.mat))
    mut_matrix_lower[i,to.keep] <- learned.sigs[to.keep,i] - sterr * qt(0.975,30)
    mut_matrix_upper[i,to.keep] <- learned.sigs[to.keep,i] + sterr * qt(0.975,30)
    mut_matrix_lower[i,to.keep][mut_matrix_lower[i,to.keep]<0] <- 0
  }
} 

learned.sigs.exome <- learned.sigs
mut_matrix_lower.exome <- t(mut_matrix_lower)
mut_matrix_upper.exome <- t(mut_matrix_upper)
for (i in 1:nrow(learned.sigs.exome)) {
  learned.sigs.exome[i,1:96] <- learned.sigs.exome[i,1:96] / trinucleotide.freq.factor.ex[types]
  mut_matrix_lower.exome[i,1:96] <- mut_matrix_lower.exome[i,1:96] + learned.sigs.exome[i,1:96] - learned.sigs[i,1:96]
  mut_matrix_upper.exome[i,1:96] <- mut_matrix_upper.exome[i,1:96] + learned.sigs.exome[i,1:96] - learned.sigs[i,1:96]
}

plot_96_profile_CI(mut_matrix = t(learned.sigs[c(2:3,5),1:96]), mut_matrix_lower = mut_matrix_lower[,c(2:3,5)], mut_matrix_upper = mut_matrix_upper[,c(2:3,5)], CI=TRUE)
plot_96_profile_CI(mut_matrix = t(learned.sigs.exome[c(2:3,5),1:96]), mut_matrix_lower = t(mut_matrix_lower.exome[c(2:3,5),]), mut_matrix_upper = t(mut_matrix_upper.exome[c(2:3,5),]), CI=TRUE)


############

random_signature <- function(sig, sig_low, sig_up) {
  res <- sapply(1:length(sig), function(i) rnorm(1,mean=sig[i],sd=(sig_up[i]-sig[i])/2))
  res[res<0] <- 0
  res[is.nan(res)] <- sig[is.nan(res)]
  return(res)
}

pdf('Histograms.pdf', 12,12)
par(mfrow=c(3,2))
hist(sapply(1:ncol(set[[5]]), function(s) cosine(learned.sigs.exome[3,1:96],set[[5]][1:96,s])), breaks=100,
     main = 'Similarity of random MMR-1 vs fixed mlh-1', xlab = 'Similarity', xlim=c(0,1))
abline(v=0.84, col='red')

hist(sapply(1:ncol(set[[5]]), function(s) cosine(learned.sigs.exome[3,-c(35,39,43,47,97:104)],set[[5]][-c(35,39,43,47,97:104),s])), breaks=100,
     main='Similarity of random MMR-1 vs fixed mlh-1 without C>T at CpG', xlab = 'Similarity', xlim=c(0,1))
abline(v=0.92, col='red')

hist(sapply(1:ncol(set[[5]]), function(s) cosine(random_signature(sig=learned.sigs.exome[3,1:96],sig_low = mut_matrix_lower.exome[3,],
                                                                  sig_up = mut_matrix_upper.exome[3,]),set[[5]][-c(97:104),s])),
     breaks=100,
     main = 'Similarity of random MMR-1 vs random mlh-1', xlab = 'Similarity', xlim=c(0,1))
abline(v=0.84, col='red')
hist(sapply(1:ncol(set[[5]]), function(s) cosine(random_signature(sig=learned.sigs.exome[3,c(1:96)[-c(35,39,43,47)]],sig_low = mut_matrix_lower.exome[3,-c(35,39,43,47)],
                                                                  sig_up = mut_matrix_upper.exome[3,-c(35,39,43,47)]),set[[5]][-c(35,39,43,47,97:104),s])),
     breaks=100,
     main = 'Similarity of random MMR-1 vs random mlh-1 without C>T at CpG', xlab = 'Similarity', xlim=c(0,1))
abline(v=0.92, col='red')
hist(sapply(1:1000, function(s) cosine(random_signature(sig=learned.sigs[2,1:96],sig_low = mut_matrix_lower[,2],sig_up = mut_matrix_upper[,2]),
                                                 random_signature(sig=learned.sigs[5,1:96],sig_low = mut_matrix_lower[,5],sig_up = mut_matrix_upper[,5]))),
     breaks=20,
     main = 'Similarity of random pole-4;pms-2 vs random mlh-1', xlab = 'Similarity', xlim=c(0,1))
hist(sapply(1:1000, function(s) cosine(random_signature(sig=learned.sigs[3,1:96],sig_low = mut_matrix_lower[,3],sig_up = mut_matrix_upper[,3]),
                                       random_signature(sig=learned.sigs[5,1:96],sig_low = mut_matrix_lower[,5],sig_up = mut_matrix_upper[,5]))),
     breaks=20,
     main = 'Similarity of random pole-4;pms-2 vs random pms-2', xlab = 'Similarity', xlim=c(0,1))
dev.off()

hist(sapply(1:1000, function(s) cosine(random_signature(sig=learned.sigs[2,1:96],sig_low = mut_matrix_lower[,2],sig_up = mut_matrix_upper[,2]),
                                       random_signature(sig=learned.sigs[3,1:96],sig_low = mut_matrix_lower[,3],sig_up = mut_matrix_upper[,3]))),
     breaks=20,
     main = 'Similarity of random pms-2 vs random mlh-1', xlab = 'Similarity', xlim=c(0,1))

###########################

# Fitting Dirichlet distribution to COSMIC similarity data
hist(cosmic[upper.tri(cosmic)], prob=T, breaks=50, main='Similarities between COSMIC signatures', xlab = 'Cosine similarity', xlim=c(0,1))
lines(density(cosmic[upper.tri(cosmic)], adjust=2), lty="dotted", col='red') 

fit.dirichlet(cancer_signatures)$p -> est.alpha2

Dr1000 <- sapply(1:30, function(x) bayesm::rdirichlet(alpha=rep(1/96,96)))
Dr1000M <- cosineM(Dr1000)
hist(Dr1000M[upper.tri(Dr1000M)],breaks=50, main='Distribution of similarities for random signatures', xlab='Similarity score', xlim=c(0,1))
abline(v = quantile(as.vector(Dr1000M[upper.tri(Dr1000M)]), 0.95), col='red', lty=2)


Dr1000 <- sapply(1:100, function(x) bayesm::rdirichlet(alpha=rep(0.05,96)))
Dr1000M <- cosineM(Dr1000)
hist(Dr1000M[upper.tri(Dr1000M)],breaks=50, main='Distribution of similarities for random signatures', xlab='Similarity score', xlim=c(0,1))
abline(v = quantile(as.vector(Dr1000M[upper.tri(Dr1000M)]), 0.95), col='red', lty=2)

sapply(1:96, function(i) var(cancer_signatures[i,]) / (mean(cancer_signatures[i])*(1 - mean(cancer_signatures))) - 1)
