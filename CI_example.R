###################################################################
### Inferring and plotting signatures with confidence intervals ###
###################################################################

source('nmSolve.R')
source('plot_profile_CI.R')
load('~/Learned_signatures.RData')
# or you can you can take your count matrix 'worm.donor.mut.mat', design matrix 'small.X' should be same
learned.sigs <- nmSolve(worm.donor.mut.mat,small.X,maxIter=10000, tol = 1e-06, div.err = 1e-10)
Y = worm.donor.mut.mat
mut_matrix_lower <- matrix(NA,nrow=ncol(learned.sigs),ncol=nrow(learned.sigs),dimnames=list(colnames(learned.sigs),row.names(learned.sigs)))
mut_matrix_upper <- matrix(NA,nrow=ncol(learned.sigs),ncol=nrow(learned.sigs),dimnames=list(colnames(learned.sigs),row.names(learned.sigs)))
library(MASS)
for (i in 1:ncol(learned.sigs)) {
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
library(ggplot2)
library(reshape2)
pdf('sigs_with_CI.pdf',15,5)
plot_96_profile_CI(mut_matrix = t(learned.sigs[c(2:3,5),1:96]), mut_matrix_lower = mut_matrix_lower[1:96,c(2:3,5)], mut_matrix_upper = mut_matrix_upper[1:96,c(2:3,5)], CI=TRUE,ymax=0.1)
dev.off()

mmr <- mut_mat[1:96,c('CD0134a','CD0134c','CD0134d','CD0135a','CD0135c','CD0135d')]
polepms <- mut_mat[1:96,c('CD0246d','CD0246e')]

plot_sig(cbind(rowMeans(mmr), rowMeans(polepms)))
low <- cbind(rowMeans(mmr) - apply(mmr,1,sd)/sqrt(6), rowMeans(polepms) - apply(polepms,1,sd)/sqrt(2) )
high <- cbind(rowMeans(mmr) + apply(mmr,1,sd)/sqrt(6) , rowMeans(polepms) + apply(polepms,1,sd)/sqrt(2) )
colnames(low) = colnames(high) = c('mmr','polepms')
ts <- cbind(rowMeans(mmr), rowMeans(polepms))
colnames(ts) <- c('mmr','polepms')
p <- plot_96_profile_CI(mut_matrix = ts,
                   mut_matrix_lower = low, mut_matrix_upper = high, CI=T, ymax=0.15)
ggsave(plot=plot, file='~/1E.pdf', width=6, height=3)
