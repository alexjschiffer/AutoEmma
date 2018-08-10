
# Supporting Function: EMMA REMLE
emma.REMLE <- function(y, X, K, Z=NULL, ngrids=100, llim=-10, ulim=10,
  esp=1e-10, eig.L = NULL, eig.R = NULL) {
  n <- length(y)
  t <- nrow(K)
  q <- ncol(X)

#  stopifnot(nrow(K) == t)
  stopifnot(ncol(K) == t)
  stopifnot(nrow(X) == n)

  if ( det(crossprod(X,X)) == 0 ) {
    warning("X is singular")
    return (list(REML=0,delta=0,ve=0,vg=0))
  }

  if ( is.null(Z) ) {
    if ( is.null(eig.R) ) {
      eig.R <- emma.eigen.R.wo.Z(K,X)
    }
    etas <- crossprod(eig.R$vectors,y)

    logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
    m <- length(logdelta)
    delta <- exp(logdelta)
    Lambdas <- matrix(eig.R$values,n-q,m) + matrix(delta,n-q,m,byrow=TRUE)
    Etasq <- matrix(etas*etas,n-q,m)
    LL <- 0.5*((n-q)*(log((n-q)/(2*pi))-1-log(colSums(Etasq/Lambdas)))-colSums(log(Lambdas)))
    dLL <- 0.5*delta*((n-q)*colSums(Etasq/(Lambdas*Lambdas))/colSums(Etasq/Lambdas)-colSums(1/Lambdas))

    optlogdelta <- vector(length=0)
    optLL <- vector(length=0)
    if ( dLL[1] < esp ) {
      optlogdelta <- append(optlogdelta, llim)
      optLL <- append(optLL, emma.delta.REML.LL.wo.Z(llim,eig.R$values,etas))
    }
    if ( dLL[m-1] > 0-esp ) {
      optlogdelta <- append(optlogdelta, ulim)
      optLL <- append(optLL, emma.delta.REML.LL.wo.Z(ulim,eig.R$values,etas))
    }

    for( i in 1:(m-1) )
      {
        if ( ( dLL[i]*dLL[i+1] < 0-esp*esp ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) )
        {
          r <- uniroot(emma.delta.REML.dLL.wo.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, etas=etas)
          optlogdelta <- append(optlogdelta, r$root)
          optLL <- append(optLL, emma.delta.REML.LL.wo.Z(r$root,eig.R$values, etas))
        }
      }
#    optdelta <- exp(optlogdelta)
  }
  else {
    if ( is.null(eig.R) ) {
      eig.R <- emma.eigen.R.w.Z(Z,K,X)
    }
    etas <- crossprod(eig.R$vectors,y)
    etas.1 <- etas[1:(t-q)]
    etas.2 <- etas[(t-q+1):(n-q)]
    etas.2.sq <- sum(etas.2*etas.2)

    logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
    m <- length(logdelta)
    delta <- exp(logdelta)
    Lambdas <- matrix(eig.R$values,t-q,m) + matrix(delta,t-q,m,byrow=TRUE)
    Etasq <- matrix(etas.1*etas.1,t-q,m)
    dLL <- 0.5*delta*((n-q)*(colSums(Etasq/(Lambdas*Lambdas))+etas.2.sq/(delta*delta))/(colSums(Etasq/Lambdas)+etas.2.sq/delta)-(colSums(1/Lambdas)+(n-t)/delta))

    optlogdelta <- vector(length=0)
    optLL <- vector(length=0)
    if ( dLL[1] < esp ) {
      optlogdelta <- append(optlogdelta, llim)
      optLL <- append(optLL, emma.delta.REML.LL.w.Z(llim,eig.R$values,etas.1,n,t,etas.2.sq))
    }
    if ( dLL[m-1] > 0-esp ) {
      optlogdelta <- append(optlogdelta, ulim)
      optLL <- append(optLL, emma.delta.REML.LL.w.Z(ulim,eig.R$values,etas.1,n,t,etas.2.sq))
    }

    for( i in 1:(m-1) )
      {
        if ( ( dLL[i]*dLL[i+1] < 0-esp*esp ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) )
        {
          r <- uniroot(emma.delta.REML.dLL.w.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, etas.1=etas.1, n=n, t1=t, etas.2.sq = etas.2.sq )
          optlogdelta <- append(optlogdelta, r$root)
          optLL <- append(optLL, emma.delta.REML.LL.w.Z(r$root,eig.R$values, etas.1, n, t, etas.2.sq ))
        }
      }
#    optdelta <- exp(optlogdelta)
  }

  maxdelta <- exp(optlogdelta[which.max(optLL)])
  maxLL <- max(optLL)
  if ( is.null(Z) ) {
    maxva <- sum(etas*etas/(eig.R$values+maxdelta))/(n-q)
  }
  else {
    maxva <- (sum(etas.1*etas.1/(eig.R$values+maxdelta))+etas.2.sq/maxdelta)/(n-q)
  }
  maxve <- maxva*maxdelta

  return (list(REML=maxLL,delta=maxdelta,ve=maxve,vg=maxva))
}

emma.eigen.L <- function(Z,K,complete=TRUE) {
  if ( is.null(Z) ) {
    return(emma.eigen.L.wo.Z(K))
  }
  else {
    return(emma.eigen.L.w.Z(Z,K,complete))
  }
}

emma.eigen.L.wo.Z <- function(K) {
  eig <- eigen(K,symmetric=TRUE)
  return(list(values=eig$values,vectors=eig$vectors))
}

emma.eigen.L.w.Z <- function(Z,K,complete=TRUE) {
  if ( complete == FALSE ) {
    vids <- colSums(Z)>0
    Z <- Z[,vids]
    K <- K[vids,vids]
  }
  eig <- eigen(K%*%crossprod(Z,Z),symmetric=FALSE,EISPACK=TRUE)
  return(list(values=eig$values,vectors=qr.Q(qr(Z%*%eig$vectors),complete=TRUE)))
}

emma.eigen.R <- function(Z,K,X,complete=TRUE) {
  if ( ncol(X) == 0 ) {
    return(emma.eigen.L(Z,K))
  }
  else if ( is.null(Z) ) {
    return(emma.eigen.R.wo.Z(K,X))
  }
  else {
    return(emma.eigen.R.w.Z(Z,K,X,complete))
  }
}

emma.eigen.R.wo.Z <- function(K, X) {
  n <- nrow(X)
  q <- ncol(X)
  S <- diag(n)-X%*%solve(crossprod(X,X))%*%t(X)
  eig <- eigen(S%*%(K+diag(1,n))%*%S,symmetric=TRUE)
  stopifnot(!is.complex(eig$values))
  return(list(values=eig$values[1:(n-q)]-1,vectors=eig$vectors[,1:(n-q)]))
}

emma.eigen.R.w.Z <- function(Z, K, X, complete = TRUE) {
  if ( complete == FALSE ) {
    vids <-  colSums(Z) > 0
    Z <- Z[,vids]
    K <- K[vids,vids]
  }
  n <- nrow(Z)
  t <- ncol(Z)
  q <- ncol(X)

  SZ <- Z - X%*%solve(crossprod(X,X))%*%crossprod(X,Z)
  eig <- eigen(K%*%crossprod(Z,SZ),symmetric=FALSE,EISPACK=TRUE)
  if ( is.complex(eig$values) ) {
    eig$values <- Re(eig$values)
    eig$vectors <- Re(eig$vectors)
  }
  qr.X <- qr.Q(qr(X))
  return(list(values=eig$values[1:(t-q)],
              vectors=qr.Q(qr(cbind(SZ%*%eig$vectors[,1:(t-q)],qr.X)),
                complete=TRUE)[,c(1:(t-q),(t+1):n)]))
}

emma.delta.ML.LL.wo.Z <- function(logdelta, lambda, etas, xi) {
  n <- length(xi)
  delta <- exp(logdelta)
  return( 0.5*(n*(log(n/(2*pi))-1-log(sum((etas*etas)/(lambda+delta))))-sum(log(xi+delta))) )
}

emma.delta.ML.LL.w.Z <- function(logdelta, lambda, etas.1, xi.1, n, etas.2.sq ) {
  t <- length(xi.1)
  delta <- exp(logdelta)
#  stopifnot(length(lambda) == length(etas.1))
  return( 0.5*(n*(log(n/(2*pi))-1-log(sum(etas.1*etas.1/(lambda+delta))+etas.2.sq/delta))-(sum(log(xi.1+delta))+(n-t)*logdelta)) )
}

emma.delta.ML.dLL.wo.Z <- function(logdelta, lambda, etas, xi) {
  n <- length(xi)
  delta <- exp(logdelta)
  etasq <- etas*etas
  ldelta <- lambda+delta
  return( 0.5*(n*sum(etasq/(ldelta*ldelta))/sum(etasq/ldelta)-sum(1/(xi+delta))) )
}

emma.delta.ML.dLL.w.Z <- function(logdelta, lambda, etas.1, xi.1, n, etas.2.sq ) {
  t <- length(xi.1)
  delta <- exp(logdelta)
  etasq <- etas.1*etas.1
  ldelta <- lambda+delta
  return( 0.5*(n*(sum(etasq/(ldelta*ldelta))+etas.2.sq/(delta*delta))/(sum(etasq/ldelta)+etas.2.sq/delta)-(sum(1/(xi.1+delta))+(n-t)/delta) ) )
}

emma.delta.REML.LL.wo.Z <- function(logdelta, lambda, etas) {
  nq <- length(etas)
  delta <-  exp(logdelta)
  return( 0.5*(nq*(log(nq/(2*pi))-1-log(sum(etas*etas/(lambda+delta))))-sum(log(lambda+delta))) )
}

emma.delta.REML.LL.w.Z <- function(logdelta, lambda, etas.1, n, t, etas.2.sq ) {
  tq <- length(etas.1)
  nq <- n - t + tq
  delta <-  exp(logdelta)
  return( 0.5*(nq*(log(nq/(2*pi))-1-log(sum(etas.1*etas.1/(lambda+delta))+etas.2.sq/delta))-(sum(log(lambda+delta))+(n-t)*logdelta)) )
}

emma.delta.REML.dLL.wo.Z <- function(logdelta, lambda, etas) {
  nq <- length(etas)
  delta <- exp(logdelta)
  etasq <- etas*etas
  ldelta <- lambda+delta
  return( 0.5*(nq*sum(etasq/(ldelta*ldelta))/sum(etasq/ldelta)-sum(1/ldelta)) )
}

emma.delta.REML.dLL.w.Z <- function(logdelta, lambda, etas.1, n, t1, etas.2.sq ) {
  t <- t1
  tq <- length(etas.1)
  nq <- n - t + tq
  delta <- exp(logdelta)
  etasq <- etas.1*etas.1
  ldelta <- lambda+delta
  return( 0.5*(nq*(sum(etasq/(ldelta*ldelta))+etas.2.sq/(delta*delta))/(sum(etasq/ldelta)+etas.2.sq/delta)-(sum(1/ldelta)+(n-t)/delta)) )
}

emma.MLE <- function(y, X, K, Z=NULL, ngrids=100, llim=-10, ulim=10,
  esp=1e-10, eig.L = NULL, eig.R = NULL)
{
  n <- length(y)
  t <- nrow(K)
  q <- ncol(X)

#  stopifnot(nrow(K) == t)
  stopifnot(ncol(K) == t)
  stopifnot(nrow(X) == n)

  if ( det(crossprod(X,X)) == 0 ) {
    warning("X is singular")
    return (list(ML=0,delta=0,ve=0,vg=0))
  }

  if ( is.null(Z) ) {
    if ( is.null(eig.L) ) {
      eig.L <- emma.eigen.L.wo.Z(K)
    }
    if ( is.null(eig.R) ) {
      eig.R <- emma.eigen.R.wo.Z(K,X)
    }
    etas <- crossprod(eig.R$vectors,y)


    logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
    m <- length(logdelta)
    delta <- exp(logdelta)
    Lambdas <- matrix(eig.R$values,n-q,m) + matrix(delta,n-q,m,byrow=TRUE)
    Xis <- matrix(eig.L$values,n,m) + matrix(delta,n,m,byrow=TRUE)
    Etasq <- matrix(etas*etas,n-q,m)
    LL <- 0.5*(n*(log(n/(2*pi))-1-log(colSums(Etasq/Lambdas)))-colSums(log(Xis)))
    dLL <- 0.5*delta*(n*colSums(Etasq/(Lambdas*Lambdas))/colSums(Etasq/Lambdas)-colSums(1/Xis))

    optlogdelta <- vector(length=0)
    optLL <- vector(length=0)
    if ( dLL[1] < esp ) {
      optlogdelta <- append(optlogdelta, llim)
      optLL <- append(optLL, emma.delta.ML.LL.wo.Z(llim,eig.R$values,etas,eig.L$values))
    }
    if ( dLL[m-1] > 0-esp ) {
      optlogdelta <- append(optlogdelta, ulim)
      optLL <- append(optLL, emma.delta.ML.LL.wo.Z(ulim,eig.R$values,etas,eig.L$values))
    }

    for( i in 1:(m-1) )
      {
        if ( ( dLL[i]*dLL[i+1] < 0-esp*esp ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) )
        {
          r <- uniroot(emma.delta.ML.dLL.wo.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, etas=etas, xi=eig.L$values)
          optlogdelta <- append(optlogdelta, r$root)
          optLL <- append(optLL, emma.delta.ML.LL.wo.Z(r$root,eig.R$values, etas, eig.L$values))
        }
      }
#    optdelta <- exp(optlogdelta)
  }
  else {
    if ( is.null(eig.L) ) {
      eig.L <- emma.eigen.L.w.Z(Z,K)
    }
    if ( is.null(eig.R) ) {
      eig.R <- emma.eigen.R.w.Z(Z,K,X)
    }
    etas <- crossprod(eig.R$vectors,y)
    etas.1 <- etas[1:(t-q)]
    etas.2 <- etas[(t-q+1):(n-q)]
    etas.2.sq <- sum(etas.2*etas.2)

    logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim

    m <- length(logdelta)
    delta <- exp(logdelta)
    Lambdas <- matrix(eig.R$values,t-q,m) + matrix(delta,t-q,m,byrow=TRUE)
    Xis <- matrix(eig.L$values,t,m) + matrix(delta,t,m,byrow=TRUE)
    Etasq <- matrix(etas.1*etas.1,t-q,m)
    #LL <- 0.5*(n*(log(n/(2*pi))-1-log(colSums(Etasq/Lambdas)+etas.2.sq/delta))-colSums(log(Xis))+(n-t)*log(deltas))
    dLL <- 0.5*delta*(n*(colSums(Etasq/(Lambdas*Lambdas))+etas.2.sq/(delta*delta))/(colSums(Etasq/Lambdas)+etas.2.sq/delta)-(colSums(1/Xis)+(n-t)/delta))

    optlogdelta <- vector(length=0)
    optLL <- vector(length=0)
    if ( dLL[1] < esp ) {
      optlogdelta <- append(optlogdelta, llim)
      optLL <- append(optLL, emma.delta.ML.LL.w.Z(llim,eig.R$values,etas.1,eig.L$values,n,etas.2.sq))
    }
    if ( dLL[m-1] > 0-esp ) {
      optlogdelta <- append(optlogdelta, ulim)
      optLL <- append(optLL, emma.delta.ML.LL.w.Z(ulim,eig.R$values,etas.1,eig.L$values,n,etas.2.sq))
    }

    for( i in 1:(m-1) )
      {
        if ( ( dLL[i]*dLL[i+1] < 0-esp*esp ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) )
        {
          r <- uniroot(emma.delta.ML.dLL.w.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, etas.1=etas.1, xi.1=eig.L$values, n=n, etas.2.sq = etas.2.sq )
          optlogdelta <- append(optlogdelta, r$root)
          optLL <- append(optLL, emma.delta.ML.LL.w.Z(r$root,eig.R$values, etas.1, eig.L$values, n, etas.2.sq ))
        }
      }
#    optdelta <- exp(optlogdelta)
  }

  maxdelta <- exp(optlogdelta[which.max(optLL)])
  maxLL <- max(optLL)
  if ( is.null(Z) ) {
    maxva <- sum(etas*etas/(eig.R$values+maxdelta))/n
  }
  else {
    maxva <- (sum(etas.1*etas.1/(eig.R$values+maxdelta))+etas.2.sq/maxdelta)/n
  }
  maxve <- maxva*maxdelta

  return (list(ML=maxLL,delta=maxdelta,ve=maxve,vg=maxva))
}

emma.MLE.noX <- function(y, K, Z=NULL, ngrids=100, llim=-10, ulim=10,
  esp=1e-10, eig.L = NULL)
{
  n <- length(y)
  t <- nrow(K)

#  stopifnot(nrow(K) == t)
  stopifnot(ncol(K) == t)

  if ( is.null(Z) ) {
    if ( is.null(eig.L) ) {
      eig.L <- emma.eigen.L.wo.Z(K)
    }
    etas <- crossprod(eig.L$vectors,y)

    logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
    m <- length(logdelta)
    delta <- exp(logdelta)
    Xis <- matrix(eig.L$values,n,m) + matrix(delta,n,m,byrow=TRUE)
    Etasq <- matrix(etas*etas,n,m)
    LL <- 0.5*(n*(log(n/(2*pi))-1-log(colSums(Etasq/Xis)))-colSums(log(Xis)))
    dLL <- 0.5*delta*(n*colSums(Etasq/(Xis*Xis))/colSums(Etasq/Xis)-colSums(1/Xis))

    optlogdelta <- vector(length=0)
    optLL <- vector(length=0)
    #print(dLL)
    if ( dLL[1] < esp ) {
      optlogdelta <- append(optlogdelta, llim)
      optLL <- append(optLL, emma.delta.ML.LL.wo.Z(llim,eig.L$values,etas,eig.L$values))
    }
    if ( dLL[m-1] > 0-esp ) {
      optlogdelta <- append(optlogdelta, ulim)
      optLL <- append(optLL, emma.delta.ML.LL.wo.Z(ulim,eig.L$values,etas,eig.L$values))
    }

    for( i in 1:(m-1) )
      {
        #if ( ( dLL[i]*dLL[i+1] < 0 ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) )
        if ( ( dLL[i]*dLL[i+1] < 0-esp*esp ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) )
        {
          r <- uniroot(emma.delta.ML.dLL.wo.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.L$values, etas=etas, xi=eig.L$values)
          optlogdelta <- append(optlogdelta, r$root)
          optLL <- append(optLL, emma.delta.ML.LL.wo.Z(r$root,eig.L$values, etas, eig.L$values))
        }
      }
#    optdelta <- exp(optlogdelta)
  }
  else {
    if ( is.null(eig.L) ) {
      eig.L <- emma.eigen.L.w.Z(Z,K)
    }
    etas <- crossprod(eig.L$vectors,y)
    etas.1 <- etas[1:t]
    etas.2 <- etas[(t+1):n]
    etas.2.sq <- sum(etas.2*etas.2)

    logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim

    m <- length(logdelta)
    delta <- exp(logdelta)
    Xis <- matrix(eig.L$values,t,m) + matrix(delta,t,m,byrow=TRUE)
    Etasq <- matrix(etas.1*etas.1,t,m)
    #LL <- 0.5*(n*(log(n/(2*pi))-1-log(colSums(Etasq/Lambdas)+etas.2.sq/delta))-colSums(log(Xis))+(n-t)*log(deltas))
    dLL <- 0.5*delta*(n*(colSums(Etasq/(Xis*Xis))+etas.2.sq/(delta*delta))/(colSums(Etasq/Xis)+etas.2.sq/delta)-(colSums(1/Xis)+(n-t)/delta))

    optlogdelta <- vector(length=0)
    optLL <- vector(length=0)
    if ( dLL[1] < esp ) {
      optlogdelta <- append(optlogdelta, llim)
      optLL <- append(optLL, emma.delta.ML.LL.w.Z(llim,eig.L$values,etas.1,eig.L$values,n,etas.2.sq))
    }
    if ( dLL[m-1] > 0-esp ) {
      optlogdelta <- append(optlogdelta, ulim)
      optLL <- append(optLL, emma.delta.ML.LL.w.Z(ulim,eig.L$values,etas.1,eig.L$values,n,etas.2.sq))
    }

    for( i in 1:(m-1) )
      {
        if ( ( dLL[i]*dLL[i+1] < 0-esp*esp ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) )
        {
          r <- uniroot(emma.delta.ML.dLL.w.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.L$values, etas.1=etas.1, xi.1=eig.L$values, n=n, etas.2.sq = etas.2.sq )
          optlogdelta <- append(optlogdelta, r$root)
          optLL <- append(optLL, emma.delta.ML.LL.w.Z(r$root,eig.L$values, etas.1, eig.L$values, n, etas.2.sq ))
        }
      }
#    optdelta <- exp(optlogdelta)
  }

  maxdelta <- exp(optlogdelta[which.max(optLL)])
  maxLL <- max(optLL)
  if ( is.null(Z) ) {
    maxva <- sum(etas*etas/(eig.L$values+maxdelta))/n
  }
  else {
    maxva <- (sum(etas.1*etas.1/(eig.L$values+maxdelta))+etas.2.sq/maxdelta)/n
  }
  maxve <- maxva*maxdelta

  return (list(ML=maxLL,delta=maxdelta,ve=maxve,vg=maxva))
}
