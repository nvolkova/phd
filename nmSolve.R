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

