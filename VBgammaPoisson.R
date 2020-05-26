#########################################################
#  Gamma-poisson NMF with variational inference         #
#  EMBL-EBI, August 2017, N.Volkova nvolkova@ebi.ac.uk  #
#########################################################

solveByNewton <- function(A, Z, N = 100, tol = 1e-02) {
  newA <- A
  for (i in 1:N) {
    delta <- (log(newA) - digamma(newA) + 1 - Z ) / (1/newA - trigamma(newA))
    oldA = newA
    negidx <- which(newA<delta,arr.ind = T)
    while(nrow(negidx)>0) {
      for (k in 1:nrow(negidx)) 
        delta[negidx[k,1],negidx[k,2]] <-  delta[negidx[k,1],negidx[k,2]] / 2
      negidx <- which(newA<delta,arr.ind = T)
    }
    newA <- newA - delta
    if (abs(mean(newA - oldA)) < tol)
      break
  }
  #if (i==N) print("Did not converge\n")
  return(newA)
}

nmSolve.vb.gamma <- function(Y, X, A, B, maxiter = 10000, tol=1e-6, div.err=1e-6, update.hyp=1) {
  
  # math. expectation of E; initialise as Gamma(A_e, B_e)
  E <- sapply(1:ncol(Y),function(i) 
    sapply(1:ncol(X), function(j) rgamma(1,shape=A[j,i],scale=B[j,i]/A[j,i])))
  # exp(ME[log(E)])
  L <- E
  # Change of parameters
  err = 0.01
  # Change of divergence
  div.change <- 2 * div.err
  
  iter=1
  # NA value mask
  mask <- !(is.na(Y))
  Y[is.na(Y)] <- 0
  # Initial divergence
  divergence.old <- mean(Y * log(Y/(X%*%E+.Machine$double.eps)) - Y + X%*%E, na.rm=T)
  
  while (iter<maxiter & (err > tol) & (abs(div.change) > div.err)) {
    
    # sufficient statistics
    sigma <- L * (t(X) %*% (Y / (X %*% E + .Machine$double.eps)))
    # Hyperparameter update
    alpha <- A + sigma
    beta <- 1 / (A/B + t(X) %*% mask)
    # change of coefficients
    err <- mean(abs(E - alpha*beta)/(E+.Machine$double.eps), na.rm=TRUE)
    # Mean update
    E <- alpha * beta
    # Compute bound:
    #bound = sum( - mask*(X %*% E) - lgamma(Y+1)) +
    #  sum(-(mask*Y) * ( ((X * log(X+.Machine$double.eps))%*%E + X%*%(E*log(E)))/(X%*%(E+.Machine$double.eps))  - log(X %*% L))) +
    #  sum(- (A/B)*E - lgamma(A) + A*log(A/B)) +
    #  sum(alpha * (log(beta) + 1) + log(alpha))
    # Exp(ME[log(E)]) update
    L <- exp(digamma(alpha)) * beta
    # Update hyperparameters:
    if (update.hyp==1) {
      Z <- E / B - log( L / B + .Machine$double.eps)
      A <- solveByNewton(A, Z)
      B <- E
    }
    divergence <- mean(Y * log(Y/(X%*%E+.Machine$double.eps)) - Y + X%*%E, na.rm=T)
    div.change <- divergence.old - divergence
    divergence.old = divergence
    
    iter = iter+1
    if(iter %% 100 == 0) cat(round(-log10(err)))
  }
  E <- alpha * beta
  divergence <- mean(Y * log(Y/(X%*%E+.Machine$double.eps)) - Y + X%*%E, na.rm=T)
  err <- mean(abs(E - alpha*beta)/(E+.Machine$double.eps), na.rm=TRUE)
  div.change <- divergence.old - divergence
  print(c(err, div.change, iter))
  return(list(alpha,beta,A,B))
  
}

KLdivergence <- function(Y, X) {
  return(mean((Y * log(Y/(X+.Machine$double.eps))) - Y + X, na.rm=T))
}

MAPnmSolve <- function(Y, X, A=0, B=1, maxIter = 10000, tol=1e-6, div.err=1e-6) {
  s = ncol(X)
  m = ncol(Y)
  if (length(which(A==0))>0) {
    E1 <- E2 <- matrix(runif(s*m,1e-3, 1), ncol = m)
  } else E1 <- E2 <- matrix(rgamma(s*m, shape = A, scale = B/A), ncol = m)
  err <- 2*tol
  mask <- !(is.na(Y))
  Y[is.na(Y)] <- 0
  
  iter <- 1
  divergence.old <-  mean(Y*log(Y/(X %*%E2 + .Machine$double.eps)) - Y + X%*%E2, na.rm=T)
  div.change <- 2 * div.err
  while (iter < maxIter & err > tol & abs(div.change) > div.err) {
    E1 <- E2
    E2 <- (A + E1 * (t(X)%*%(Y/(X%*%E1 +.Machine$double.eps)))) / (A/B + t(X)%*%mask +.Machine$double.eps)
    iter <- iter + 1
    err <- mean(abs(E2 - E1)/(E1+.Machine$double.eps), na.rm=TRUE)
    divergence <- mean(Y*log(Y/(X %*%E2 + .Machine$double.eps)) - Y + X%*%E2, na.rm=T) # KL distance from D to P%*%E2
    div.change <- divergence.old - divergence
    divergence.old = divergence
    if(iter %% 100 == 0) cat(round(-log10(err)))
    # add likelihood convergence
  }
  cat("\n")
  if(iter == maxIter) warning(paste("No convergence after",iter, "iterations."))
  print(c(err, div.change, divergence, iter))
  return(E2)
}

################################################################################

# MCMC sampling

Gibbs_sampler <- function(Y, X, A=NA, B=NA, maxiter = 1000, tol=1e-6, div.err=1e-6) {
  
  listE <- list()
  
  if(is.na(A))
    A <- matrix(0.01,nrow(E),ncol(E))
  if(is.na(B))
    B <- matrix(0.01,nrow(E),ncol(E))
  # Initialize E as Gamma RV from A and B
  E <- sapply(1:ncol(Y),function(i) 
    sapply(1:ncol(X), function(j) rgamma(1,shape=A[j,i],rate=B[j,i])))
  # Change of parameters
  err = 100
  # Change of divergence
  div.change <- 2 * div.err
  # Iterations
  iter=1
  # NA value mask
  mask <- !(is.na(Y))
  Y[is.na(Y)] <- 0
  # Initial divergence
  divergence.old <- mean(Y * log(Y/(X%*%E+.Machine$double.eps)) - Y + X%*%E, na.rm=T)

  plot(NA, NA, xlim=c(1,maxiter), ylim=c(-20,20))
  abline(h = log10(ref.err), lty=2, col='black')
  abline(h = log10(ref.div), lty=2, col='red')
  
  while ( iter < maxiter & 
          err > tol & 
          abs(div.change) > div.err) { 
    Sys.sleep(.01)
    # Sample sources
    P <- array(sapply(1:nrow(Y),function(i) 
           sapply(1:ncol(Y), function(j)
             X[i,] * E[,j] / sum( X[i,] * E[,j] ))), dim=c(nrow(E),ncol(Y),nrow(Y)))
    # Draw latent factors
    S <- array(sapply(1:ncol(Y), function(j) {
      sapply(1:nrow(Y), function(i) rmultinom(n=ncol(X), size = Y[i,j], prob = P[,j,i]))
    }), dim=c(nrow(E),nrow(Y),ncol(Y)))
    # Sufficient statistics
    Sigma <- apply(S,c(1,3),sum)
    # Sample effects
    Anew <- A + Sigma
    Bnew <- B + t(X) %*% mask
    E <- sapply(1:ncol(Y),function(i) 
      sapply(1:ncol(X), function(j) rgamma(1,shape=Anew[j,i],rate=Bnew[j,i])))
    # Calculate divergence and error
    err <- sum((Y - X%*%E)**2)
    points(iter,log10(err),pch=16,cex=0.5)
    divergence <- mean(Y * log(Y/(X%*%E+.Machine$double.eps)) - Y + X%*%E, na.rm=T)
    div.change <- divergence.old - divergence
    divergence.old <- divergence
    points(iter,log10(divergence.old),pch=16,col='red',cex=0.5)
    #points(iter,divergence.old,pch=16,cex=0.5)
    
    iter = iter+1
    if(iter %% 10 == 0) {
      cat(round(-log10(err)))
      cat(' ')
      cat(round(div.change,2))
      cat('\n')
    }
    
    if (iter>500) listE[[iter-500]] <- E
  }
  
  print(c(err, div.change, iter))
  return(list(E,alpha,beta))
  
}

Gibbs_sampler_Zhou <- function(Y, 
                               X, 
                               R = matrix(1,nrow=nrow(E),ncol=ncol(E)), 
                               P = matrix(0.1,nrow=nrow(E),ncol=ncol(E)),
                               maxiter = 1000, tol=1e-6, div.err=1e-6,
                               eps = 0.01, c = 1, c0 = 1, r0 = 1) {
  
  N = ncol(Y)
  # Initialize E as Gamma RV from A and B
  E <- t(sapply(1:ncol(X), function(j) 
    sapply(1:N, function(i) rgamma(1,shape=R[j,i],scale=P[j,i] / (1-P[j,i])))))
  # Change of parameters
  err = 100
  # Change of divergence
  div.change <- 2 * div.err
  # Iterations
  iter=1
  # NA value mask
  mask <- !(is.na(Y))
  Y[is.na(Y)] <- 0
  # Initial divergence
  divergence.old <- mean(Y * log(Y/(X%*%E+.Machine$double.eps)) - Y + X%*%E, na.rm=T)
  #par(mfrow=c(2,1))
  plot(NA, NA, xlim=c(1,maxiter), ylim=c(-20,20))
  points(1,log10(err),pch=16,cex=0.5)
  #plot(NA, NA, xlim=c(1,maxiter), ylim=c(0,2*divergence.old))
  #points(1,divergence.old,pch=16,cex=0.5)
  
  while ( iter < maxiter & 
          err > tol & 
          abs(div.change) > div.err) { 
    
    # For plotting
    Sys.sleep(.01)
    
    # Sample sources
    Pmult <- array(sapply(1:nrow(Y),function(i) 
      sapply(1:ncol(Y), function(j)
        X[i,] * E[,j] / sum( X[i,] * E[,j] ))), dim=c(nrow(E),ncol(Y),nrow(Y)))
    
    # Draw latent factors
    S <- array(sapply(1:ncol(Y), function(j) {
      sapply(1:nrow(Y), function(i) rmultinom(n=ncol(X), size = Y[i,j], prob = Pmult[,j,i]))
    }), dim=c(nrow(E),nrow(Y),ncol(Y)))
    
    # Sufficient statistics
    #Sigma <- apply(S,c(1,3),sum)
    
    # Sample effects
    Pnew <- sapply(1:nrow(E), function(k) 
      sapply(1:ncol(E), function(i) rbeta(1,c*eps + sum(S[k,,i]),c*(1-eps)+N*R[k,i])))
    R <- sapply(1:nrow(E), function(k) 
      sapply(1:ncol(E), function(i) rgamma(1, shape=c0*r0, scale=1/(c0-N*log(1-P[k,i])))))
    P <- Pnew
    E <- sapply(1:ncol(Y),function(i) 
      sapply(1:ncol(X), function(j) rgamma(1,shape=R[j]+sum(S[j,i,]),scale=P[j])))
    
    # Calculate divergence and error
    err <- sum((Y - X%*%E)**2)
    points(iter,log10(err),pch=16,cex=0.5)
    divergence <- mean(Y * log(Y/(X%*%E+.Machine$double.eps)) - Y + X%*%E, na.rm=T)
    div.change <- divergence.old - divergence
    divergence.old <- divergence
    points(iter,log10(divergence.old),pch=16,col='red',cex=0.5)
    
    
    iter = iter+1
    if(iter %% 10 == 0) {
      cat(round(-log10(err)))
      cat(' ')
      cat(round(div.change,2))
      cat('\n')
    }
  }
  
  print(c(err, div.change, iter))
  return(list(E,alpha,beta))
  
}
# Estimate the marginal likelihood P(Y|A,B) to be able to update A,B?
#Chibbs_marginal_likelihood <- function() {
#  
#}

