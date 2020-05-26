# estimate per-experiment variance and mean (quiasi-poisson parameter)

table(substr(worms,1,6))
proper.exp <- names(which(table(substr(worms,1,6))>2 & table(substr(worms,1,6))<5))

alphas <- list()
for (i in 1:length(proper.exp)) {
  test <- worms[grep(proper.exp[i],worms)]
  if (length(grep("b",test))!=0)
    test <- test[-grep("b",test)]
  
  #alphas[[i]] <- sapply(1:10, function(j) sum((spectrum[test,j]-mean(Fitted[test,j]))**2,na.rm=T) / (mean(Fitted[test,j]) + .Machine$double.eps) / length(test))
  
  lambda <- colMeans(Y[test,], na.rm=T)
  alphas[[i]] <- sapply(1:10, function(j) {
    if (is.na(lambda[j])) return(1)
    if (lambda[j] < 0.1) return(1)
    y <- sum((Y[test,j] - lambda[j])**2, na.rm=T) / (lambda[j] + .Machine$double.eps) / (length(test)-1)
    if (y>=1) return(y)
    else return(1)
  })
  
  
}
par.list <- do.call("rbind",alphas)
hist(par.list, breaks=1000, xlim=c(0,20))




qqplot(x = do.call("c",sapply(proper.exp, function(x) {
    test <- worms[grep(x,worms)] 
    if (length(grep("b",test))!=0)
      test <- test[-grep("b",test)]
  
    lambda <- colMeans(Y[test,], na.rm=T)
    alphas <- sapply(1:10, function(j) {
      if (is.na(lambda[j])) return(1)
      if (lambda[j] < 0.1) return(1)
      y <- sum((Y[test,j] - lambda[j])**2, na.rm=T) / (lambda[j] + .Machine$double.eps) / (length(test)-1)
      if (y>=1) return(y)
      else return(1)
    })
    randoms <- sapply(1:10, function(i) {
      if (alphas[i]==1)
        rpois(n=length(test), lambda[i]) - Fitted[test,i]
      else
        rqpois(n=length(test), mu=lambda[i], theta=alphas[i]) - Fitted[test,i]
    })
    randoms[is.na(randoms)] <- 0
    return(as.vector(randoms))
  })),
  y = do.call("c",sapply(proper.exp, function(x) { 
    test <- worms[grep(x,worms)]
    if (length(grep("b",test))==0)
      test <- test[-grep("b",test)]
    return(as.vector(Y[test,] - Fitted[test,]))
  })),
  main = "QQ-plot of random QPoisson quantiles vs real residuals",
  ylab = "Theoretical quantiles", xlab = "Model residuals")
abline(a=0,b=1,col="red")

qqplot(x = do.call("c",sapply(proper.exp, function(x) {
  test <- worms[grep(x,worms)]
  if (length(grep("b",test))!=0)
    test <- test[-grep("b",test)]
  lambda <- colMeans(Fitted[test,])
  randoms <- sapply(1:10, function(j) rnbinom(n = length(test), mu = lambda[j], size = lambda[j]/10 + 1) - lambda[j])
  return(as.vector(randoms))
})),
y = do.call("c",sapply(proper.exp, function(x) { 
  test <- worms[grep(x,worms)]
  if (length(grep("b",test))==0)
    test <- test[-grep("b",test)]
  return(as.vector(Y[test,] - Fitted[test,]))
})),
main = "QQ-plot of random QPoisson quantiles vs real residuals",
ylab = "Theoretical quantiles", xlab = "Model residuals")
abline(a=0,b=1,col="red")


plot(x = as.vector(Fitted), y = as.vector(Y), cex = 0.4, main="Predicted values vs observed \n all interaction", xlab = "Predicted", ylab = "Observed")
abline(a=0,b=1,col="red",lty=2)
#ind <- which(abs(Fitted - Y)>50, arr.ind = T)
#text(x = Fitted[ind], y = Y[ind], labels=CD2Mutant[rownames(ind)], pos=1, cex=0.6)

for (i in 1:10) {
  if (length(underestimated.outliers[[i]])>0) {
    points(x = Fitted[underestimated.outliers[[i]],i], y = Y[underestimated.outliers[[i]],i], col="red",pch=16)
    text(x = Fitted[underestimated.outliers[[i]],i], y = Y[underestimated.outliers[[i]],i], labels=CD2Mutant[underestimated.outliers[[i]]], pos=1, cex=0.6)
  }
  if (length(overestimated.outliers[[i]])>0) {
    points(x = Fitted[overestimated.outliers[[i]],i], y = Y[overestimated.outliers[[i]],i], col="red",pch=16)
    text(x = Fitted[overestimated.outliers[[i]],i], y = Y[overestimated.outliers[[i]],i], labels=CD2Mutant[overestimated.outliers[[i]]], pos=1, cex=0.6)
  }
}

plot(x = as.vector(Y), y = as.vector(Y) - as.vector(Fitted), cex = 0.4, main="Residuals vs observed values \n all interaction", xlab = "Observed", ylab = "Residuals")
abline(h=0,col="red",lty=2)
abline(a=0,b=1,col="red",lty=2)
abline(a=0,b=-1,col="red",lty=2)
for (i in 1:10) {
  if (length(underestimated.outliers[[i]])>0) {
    points(x = Y[underestimated.outliers[[i]],i], y = Y[underestimated.outliers[[i]],i]-Fitted[underestimated.outliers[[i]],i], col="red",pch=16)
    text(x = Y[underestimated.outliers[[i]],i], y = Y[underestimated.outliers[[i]],i]-Fitted[underestimated.outliers[[i]],i], labels=CD2Mutant[underestimated.outliers[[i]]], pos=1, cex=0.6)
  }
  if (length(overestimated.outliers[[i]])>0) {
    points(x = Y[overestimated.outliers[[i]],i], y = Y[overestimated.outliers[[i]],i]-Fitted[overestimated.outliers[[i]],i], col="red",pch=16)
    text(x = Y[overestimated.outliers[[i]],i], y = Y[overestimated.outliers[[i]],i]-Fitted[overestimated.outliers[[i]],i], labels=CD2Mutant[overestimated.outliers[[i]]], pos=1, cex=0.6)
  }
}


plot(x = as.vector(Fitted), y = as.vector(Y) - as.vector(Fitted), cex = 0.4, main="Residuals vs predicted values \n all interaction", xlab = "Predicted", ylab = "Residuals")
abline(h=0,col="red",lty=2)
abline(a=0,b=1,col="red",lty=2)
abline(a=0,b=-1,col="red",lty=2)
for (i in 1:10) {
  if (length(underestimated.outliers[[i]])>0) {
    points(x = Fitted[underestimated.outliers[[i]],i], y = Y[underestimated.outliers[[i]],i]-Fitted[underestimated.outliers[[i]],i], col="red",pch=16)
    text(x = Fitted[underestimated.outliers[[i]],i], y = Y[underestimated.outliers[[i]],i]-Fitted[underestimated.outliers[[i]],i], labels=CD2Mutant[underestimated.outliers[[i]]], pos=1, cex=0.6)
  }
  if (length(overestimated.outliers[[i]])>0) {
    points(x = Fitted[overestimated.outliers[[i]],i], y = Y[overestimated.outliers[[i]],i]-Fitted[overestimated.outliers[[i]],i], col="red",pch=16)
    text(x = Fitted[overestimated.outliers[[i]],i], y = Y[overestimated.outliers[[i]],i]-Fitted[overestimated.outliers[[i]],i], labels=CD2Mutant[overestimated.outliers[[i]]], pos=1, cex=0.6)
  }
}


group <- list()
for (type in unique(CD2Mutant)) {
  group[[type]] <- names(which(CD2Mutant==type))
}

barplot(sapply(group,length), xaxt="n", xlab="Group")
axis(side=1, las=2, at=1:length(group),labels=names(group))

interesting <- group[which(sapply(group,length)>3)]
alphas <- list()
alphas2 <- list()
for (i in 1:length(interesting)) {

  lambda <- colMeans(spectrum[interesting[[i]],], na.rm=T)
  alphas[[i]] <- sapply(1:10, function(j) sum((spectrum[interesting[[i]],j] - lambda[j])**2, na.rm=T) / lambda[j] / length(interesting[[i]]))
  alphas2[[i]] <- sapply(1:10, function(j) sum((spectrum[interesting[[i]],j] - mean(Fitted[interesting[[i]],j]))**2, na.rm=T) / mean(Fitted[interesting[[i]],j]) / length(interesting[[i]]))
    
}
par.list <- do.call("rbind",alphas)
par.list2 <- do.call("rbind",alphas2)
hist(par.list, breaks=1000, main="Index of dispersion per groups of samples")
hist(par.list2, breaks=10000,xlim=c(0,50))





table(substr(worms,1,6))
proper.exp <- names(which(table(substr(worms,1,6))>3 & table(substr(worms,1,6))<5))
alphas <- list()
alphas2 <- list()
for (i in 1:length(proper.exp)) {
  test <- worms[grep(proper.exp[i],worms)]
  test <- test[-grep("b",test)]
  
  lambda <- colMeans(spectrum[test,], na.rm=T)
  alphas[[i]] <- sapply(1:10, function(j) sum((spectrum[test,j] - lambda[j])**2, na.rm=T) / lambda[j] / (length(test)-1))
  alphas2[[i]] <- sapply(1:10, function(j) sum((spectrum[test,j] - mean(Fitted[test,j]))**2, na.rm=T) / mean(Fitted[test,j]) / (length(test)-1))
  
}
par.list <- do.call("rbind",alphas)
par.list2 <- do.call("rbind",alphas2)
par(mfrow=c(1,2))
hist(par.list, breaks=1000, main="Index of dispersion per experiment", xlab="Variance / mean")
hist(par.list2, breaks=1000, main="Index of dispersion per experiment for predicted mean", xlab="Variance / mean")


for (i in 1:10) {
  
  cov.mat <- poisI(X,E.wtgm.int[,i],Y[,i])
  to.keep <- colnames(cov.mat)[which(diag(cov.mat)!=0)]
  if (length(to.keep)!=0) {
    cov.mat <- cov.mat[to.keep,to.keep]
    SVD <- svd(cov.mat)$d
    if (length(which(SVD<1e-6))>0)
      cov.mat <- ginv(cov.mat)
    if (length(which(SVD<1e-6))==0)
      cov.mat <- solve(cov.mat)
    upper <- Fitted[,i] + qnorm(.975) * sqrt(diag(X[,to.keep] %*% cov.mat %*% t(X[,to.keep])))
    lower <- Fitted[,i] - qnorm(.975) * sqrt(diag(X[,to.keep] %*% cov.mat %*% t(X[,to.keep])))
  } else {
    upper <- NULL
    lower <- NULL
  }

}


alpha = matrix(1, nrow = nrow(spectrum), ncol = ncol(spectrum), dimnames = dimnames(spectrum))
for (name in rownames(spectrum)) {
  
  test <- names(CD2Mutant)[CD2Mutant==CD2Mutant[name]]
  if (length(test)>2)
    alpha[name,] <- sapply(1:10, function(j) {
      if (is.na(mean(spectrum[test,j],na.rm=T))) return(1)
      if (mean(spectrum[test,j],na.rm=T)>0)
        return(sum((spectrum[test,j] - mean(spectrum[test,j],na.rm=T))**2, na.rm=T) / mean(spectrum[test,j],na.rm=T) / (length(test)-1))
      else return(1)
    })
}

Y <- as.matrix(cur.spectrum)
X <- as.matrix(cur.design)
E <- nmSolveQP(D = Y, P = X, alpha=alpha)

########################################

alpha <- matrix(1, nrow=nrow(spectrum), ncol=ncol(spectrum), dimnames = dimnames(spectrum))
for (name in worms) {
  test <- worms[grep(substr(name,1,6),worms)]
  if (length(grep("b",test))!=0)
    test <- test[-grep("b",test)]
  
  if (length(test)<3)
    next
  
  lambda <- colMeans(spectrum[test,], na.rm=T)
  alpha[name,] <- sapply(1:10, function(j) {
    if (is.na(lambda[j])) return(1)
    if (lambda[j] < 0.1) return(1)
    y <- sum((spectrum[test,j] - lambda[j])**2, na.rm=T) / (lambda[j] + .Machine$double.eps) / (length(test)-1)
    if (y>=1) return(y)
    else return(1)
  })
  
}

Z <- nmSolveQP(D = as.matrix(spectrum[row.names(X),]), P = X, alpha = alpha)
Fitted.alpha <- X %*% Z

plot(as.vector(Fitted.alpha), as.vector(Fitted), main="Pois vs QPois")
abline(a=0,b=1,col="red")

plot(as.vector(Y), as.vector(Y) - as.vector(Fitted.alpha))
plot(as.vector(Y), as.vector(Y) - as.vector(Fitted))

plot(as.vector(as.matrix(spectrum)), as.vector(Fitted.alpha)*alpha)
abline(a=0,b=1,col="red")

plot(as.vector(as.matrix(spectrum)), as.vector(Fitted))
abline(a=0,b=1,col="red")

alpha <- matrix(1, nrow=nrow(spectrum), ncol=ncol(spectrum), dimnames = dimnames(spectrum))
for (name in worms) {
  test <- worms[grep(substr(name,1,6),worms)]
  if (length(grep("b",test))!=0)
    test <- test[-grep("b",test)]
  
  if (length(test)<3)
    next
  
  lambda <- colMeans(Fitted.alpha[test,], na.rm=T)
  alpha[name,] <- sapply(1:10, function(j) {
    if (is.na(lambda[j])) return(1)
    if (lambda[j] < 0.1) return(1)
    y <- sum((spectrum[test,j] - lambda[j])**2, na.rm=T) / (lambda[j] + .Machine$double.eps) / (length(test)-1)
    if (y>=1) return(y)
    else return(1)
  })
  
}
hist(alpha, breaks=10000, xlim=c(0,20))


##################################################

sigs <- t(matrix(runif(30,0,0.3),nrow=10))
for (i in 1:3)
  sigs[i,] <- sigs[i,] / sum(sigs[i,])
means <- runif(1000,10,500)

mu <- matrix(NA,nrow=1000,ncol=10)
sim_pois <- matrix(NA,nrow=1000,ncol=10)
sim_od <- matrix(NA,nrow=1000,ncol=10)
exp <- matrix(NA, nrow=1000,ncol=3)


alpha <- rbinom(n = 1000, prob = 0.8, size=1)
alpha[alpha==0] <- runif(n=length(which(alpha==0)), 3,10)
for (i in 1:1000) {
  contrib <- runif(3)
  contrib <- contrib / sum(contrib)
  exp[i,] <- contrib * means[i]
  mu[i,] <- (contrib * means[i]) %*% sigs
  sim_pois[i,] <- rpois(n=10, lambda=mu[i,])
  if (alpha[i]>1)
    sim_od[i,] <- rnbinom(n=10, mu=mu[i,], size=mu[i,]/(alpha[i]-1))
  else sim_od[i,] <- rpois(n=10, lambda=mu[i,])
}

alphas <- matrix(rep(alpha,10), nrow=1000)

###########################################
## No noise at all
restore_no_noise <- nmSolve(D = mu, P = exp)
plot(as.vector(mu), as.vector(exp %*% restore_no_noise), main = "No noise", xlab = "Observed", ylab="Predicted")
abline(a=0,b=1,col="red",lty=2)

sum((sigs-restore)**2) # 0.046
plot(as.vector(mu), as.vector(mu) - as.vector(exp %*% restore_no_noise))
###########################################
## Only Poisson noise
restore_pois <- nmSolve(D = sim_pois, P = exp)
plot(as.vector(sim_pois), as.vector(exp %*% restore_pois), main = "Poisson noise", xlab = "Observed", ylab="Predicted")
abline(a=0,b=1,col="red",lty=2)

sum((sigs-restore_pois)**2) # 0.0002
## Residuals
plot(as.vector(sim_pois), as.vector(sim_pois) - as.vector(exp %*% restore_pois))
abline(a=0, b=1, col="red")
abline(a=0,b=-1,col="red")

## Likelihood?
x <- ppois(sim_pois, exp %*% restore_pois)
x[which(exp %*% restore_pois==0)] <- 0
nx <- qnorm(x)
if (length(which(nx==Inf))>0) {
  names(which(nx==Inf)) -> infinite.tmp
  nx <- nx[-which(nx==Inf)]
}
hist(nx,breaks=100)
mean(nx)
sd(nx)
skewness(as.vector(nx))
###########################################
## Plus overdispersion
restore_overdispersed <- nmSolve(D = sim_od, P = exp)
plot(as.vector(sim_od), as.vector(exp %*% restore_overdispersed), main = "Pois noise + overdispersion", xlab = "Observed", ylab="Predicted")
abline(a=0,b=1,col="red",lty=2)
points(as.vector(sim_od[alpha>1,]), as.vector((exp %*% restore_overdispersed)[alpha>1,]), pch=18, col="blue")
points(as.vector(sim_od[alpha==1,]), as.vector((exp %*% restore_overdispersed)[alpha==1,]), col="green")

sum((sigs-restore_overdispersed)**2) # 0.0005
## Residuals
plot(as.vector(sim_od), as.vector(sim_od) - as.vector(exp %*% restore_overdispersed))
points(as.vector(sim_od[alpha>1,]), as.vector(sim_od[alpha>1,]) - as.vector((exp %*% restore_overdispersed)[alpha>1,]), pch=18, col="blue")
points(as.vector(sim_od[alpha==1,]), as.vector(sim_od[alpha==1,]) - as.vector((exp %*% restore_overdispersed)[alpha==1,]), col="green")
abline(a=0, b=1, col="red")
abline(a=0,b=-1,col="red")

## Likelihood?
x <- ppois(sim_od, exp %*% restore_overdispersed)
nx <- qnorm(x)
if (length(which(nx==Inf))>0) {
  names(which(nx==Inf)) -> infinite.tmp
  nx <- nx[-which(nx==Inf)]
}
h <- hist(nx,breaks=100)
plot(h)
lines(hist(rnorm(length(nx), mean = , sd=))$counts)

mean(nx)
sd(nx)
skewness(nx)
###########################################

# Simulate replicates with the same exposures

sigs <- t(matrix(runif(30,0,0.5),nrow=10))
for (i in 1:3)
  sigs[i,] <- sigs[i,] / sum(sigs[i,])

exp <- matrix(NA, nrow=33,ncol=3)
exp[1,] <- c(1,0,0)
for (i in 2:11)
  exp[i,] <- c(20,0,0)
exp[12,] <- c(0,1,0)
for (i in 13:22)
  exp[i,] <- c(0,20,0)
exp[23,] <- c(1,1,1)
for (i in 24:33)
  exp[i,] <- c(20,20,20)
mu <- exp %*% sigs

sim_pois <- matrix(NA,nrow=33,ncol=10)
sim_od <- matrix(NA,nrow=33,ncol=10)

alpha <- c(1,rep(1,5),rep(5,5),1,rep(1,5),rep(5,5),1,rep(1,5),rep(5,5))
for (i in 1:33) {
  sim_pois[i,] <- rpois(n=10, lambda=mu[i,])
  if (alpha[i]>1)
    sim_od[i,] <- rnbinom(n=10, mu=mu[i,], size=mu[i,]/(alpha[i]-1))
  else sim_od[i,] <- rpois(n=10, lambda=mu[i,])
}

alphas <- matrix(rep(alpha,10), nrow=33)


restore_no_noise <- nmSolve(D = mu, P = exp)
plot(as.vector(mu), as.vector(exp %*% restore_no_noise), main = "No noise", xlab = "Observed", ylab="Predicted")
abline(a=0,b=1,col="red",lty=2)
sum((sigs-restore)**2) # 0.12
plot(as.vector(mu), as.vector(mu) - as.vector(exp %*% restore_no_noise))

restore_pois <- nmSolve(D = sim_pois, P = exp)
plot(as.vector(sim_pois), as.vector(exp %*% restore_pois), main = "Poisson noise", xlab = "Observed", ylab="Predicted")
abline(a=0,b=1,col="red",lty=2)
sum((sigs-restore_pois)**2) # 0.03
plot(as.vector(sim_pois), as.vector(sim_pois) - as.vector(exp %*% restore_pois))

restore_overdispersed <- nmSolve(D = sim_od, P = exp)
plot(as.vector(sim_od), as.vector(exp %*% restore_overdispersed), main = "Pois noise + overdispersion", xlab = "Observed", ylab="Predicted")
abline(a=0,b=1,col="red",lty=2)
points(as.vector(sim_od[alpha>1,]), as.vector((exp %*% restore_overdispersed)[alpha>1,]), pch=18, col="blue")
points(as.vector(sim_od[alpha==1,]), as.vector((exp %*% restore_overdispersed)[alpha==1,]), col="green")
sum((sigs-restore_overdispersed)**2) # 0.05
plot(as.vector(sim_od), as.vector(sim_od) - as.vector(exp %*% restore_overdispersed))
points(as.vector(sim_od[alpha>1,]), as.vector(sim_od[alpha>1,]) - as.vector((exp %*% restore_overdispersed)[alpha>1,]), pch=18, col="blue")
points(as.vector(sim_od[alpha==1,]), as.vector(sim_od[alpha==1,]) - as.vector((exp %*% restore_overdispersed)[alpha==1,]), col="green")

