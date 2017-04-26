nmSolve <- function(D, P, maxIter = 10000, tol=1e-5, div.err=1e-7) {
  n <- nrow(D)
  m <- ncol(D)
  s <- ncol(P)
  tP <- t(P)
  rP <- rep(colSums(P), m)
  D <- as.matrix(D)
  P <- as.matrix(P)
  E1 <- E2 <- matrix(runif(s * m, 1e-3, 1), ncol = m)
  err <- 2*tol
  
  iter <- 1
  divergence.old <- mean(D * log(D/(P%*%E2)) - D + P%*%E2, na.rm=T)
  div.change <- 2 * div.err
  while (iter < maxIter & err > tol & abs(div.change) > div.err) {
    E1 <- E2
    E2 <- E1 * (tP %*% (D/(P %*% (E1 + .Machine$double.eps))))/rP
    iter <- iter + 1
    err <- mean(abs(E2 - E1)/(E1+.Machine$double.eps), na.rm=TRUE)
    divergence <- mean(D*log(D/(P %*% (E2 + .Machine$double.eps))) - D + P%*%E2, na.rm=T) # KL distance from D to P%*%E2
    div.change <- divergence.old - divergence
    divergence.old = divergence
    if(iter %% 100 == 0) cat(round(-log10(err)))
    # add likelihood convergence
  }
  cat("\n")
  if(iter == maxIter) warning(paste("No convergence after",iter, "iterations."))
  E2
}

nmFit <- function(D, E, maxIter = 10000, tol=1e-5, div.err=1e-7) {
  n <- nrow(D)
  m <- ncol(D)
  s <- nrow(E)
  tE <- t(E)
  rE <- rep(rowSums(E),each=n)
  D <- as.matrix(D)
  E <- as.matrix(E)
  P1 <- P2 <- matrix(runif(s * n, 1e-3, 1), nrow = n)
  err <- 2*tol
  
  iter <- 1
  divergence.old <- mean(D * log(D/(P2%*%(E + .Machine$double.eps))) - D + P2%*%E, na.rm=T)
  div.change <- 2 * div.err
  while (iter < maxIter & err > tol & abs(div.change) > div.err) {
    P1 <- P2
    P2 <- P1 * ((D/((P1 +.Machine$double.eps) %*% (E +.Machine$double.eps))) %*% tE) / (rE+.Machine$double.eps)
    iter <- iter + 1
    err <- mean(abs(P2 - P1)/(P1+.Machine$double.eps), na.rm=TRUE)
    divergence <- mean(D*log(D/((P2 + .Machine$double.eps)%*%(E+.Machine$double.eps))) - D + P2%*%E, na.rm=T) # KL distance from D to P%*%E
    div.change <- divergence.old - divergence
    divergence.old = divergence
    if(iter %% 100 == 0) cat(round(-log10(err)))
    # add likelihood convergence
  }
  cat("\n")
  if(iter == maxIter) warning(paste("No convergence after",iter, "iterations."))
  P2
}

ssnmf <- function(D, S, s=0, maxIter=100, minE = 0, whichS = 1:ncol(D)){
    n <- nrow(D)
    o <- ncol(S)
    m <- ncol(D)
    P <- cbind(S, if(s>0) matrix(runif(n*s,0,1), ncol=s) else NULL)
    P <- P/rep(colSums(P), each=nrow(P))
    colnames(P) <- c(colnames(S), if(s >0) paste0("N.",1:s) else NULL)
    E <- matrix(runif((s+o)*m, 0,1), ncol=m)
    E <- E * (t(P)%*% (D / (P %*% E))) / rep(colSums(P), m)
    
    iter <- 1
    while(iter < maxIter){
        P <- P * ((D / (P %*% E)) %*% t(E)) / rep(rowSums(E), each=n)
        P[,1:o] <- S
        P <- P/rep(colSums(P), each=nrow(P))
        E <- E * (t(P)%*% (D / (P %*% E))) / rep(colSums(P), m)
        E[E/rep(colSums(E), each=nrow(E)) < minE] <- 0
        E[-(1:o),setdiff(1:ncol(E), whichS)] <- 0
        E <- E * rep(colSums(D)/colSums(E), each=nrow(E))
        iter <- iter + 1
    }
    list(E=E,P=P)
}
