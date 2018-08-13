# EMMA using ML
emma.ML.LRT <- function(ys, xs, K, Z=NULL, X0 = NULL, ngrids=100, llim=-10, ulim=10, esp=1e-10, ponly = FALSE) {
  if ( is.null(dim(ys)) || ncol(ys) == 1 ) {
    ys <- matrix(ys,1,length(ys))
  }
  if ( is.null(dim(xs)) || ncol(xs) == 1 ) {
    xs <- matrix(xs,1,length(xs))
  }
  if ( is.null(X0) ) {
    X0 <- matrix(1,ncol(ys),1)
  }

  g <- nrow(ys)
  n <- ncol(ys)
  m <- nrow(xs)
  t <- ncol(xs)
  q0 <- ncol(X0)
  q1 <- q0 + 1

  if ( !ponly ) {
    ML1s <- matrix(nrow=m,ncol=g)
    ML0s <- matrix(nrow=m,ncol=g)
    vgs <- matrix(nrow=m,ncol=g)
    ves <- matrix(nrow=m,ncol=g)
  }
  stats <- matrix(nrow=m,ncol=g)
  ps <- matrix(nrow=m,ncol=g)
  ML0 <- vector(length=g)

  stopifnot(nrow(K) == t)
  stopifnot(ncol(K) == t)
  stopifnot(nrow(X0) == n)

  if ( sum(is.na(ys)) == 0 ) {
    eig.L <- emma.eigen.L(Z,K)
    eig.R0 <- emma.eigen.R(Z,K,X0)

    # Loops through each of the response variables (phenotypes)
    for(i in 1:g) {
      ML0[i] <- emma.MLE(ys[i,],X0,K,Z,ngrids,llim,ulim,esp,eig.L,eig.R0)$ML
    }

    x.prev <- vector(length=0)

    # Initialize Progress Bar
    pb <- txtProgressBar(min = 0, initial = 0, max = m, style = 3)

    # Loops through each of the SNPs
    for(i in 1:m) {
      vids <- !is.na(xs[i,])
      nv <- sum(vids)
      xv <- xs[i,vids]

      # Update Progress Bar
      if (i %% 10 == 0) {
        Sys.sleep(0.001)
        setTxtProgressBar(pb, i)
      }

      if ( ( mean(xv) <= 0 ) || ( mean(xv) >= 1 ) ) {
        if (!ponly) {
          stats[i,] <- rep(NA,g)
          vgs[i,] <- rep(NA,g)
          ves[i,] <- rep(NA,g)
          ML1s[i,] <- rep(NA,g)
          ML0s[i,] <- rep(NA,g)
        }
        ps[i,] = rep(1,g)
      }
      else if ( identical(x.prev, xv) ) {
        if ( !ponly ) {
          stats[i,] <- stats[i-1,]
          vgs[i,] <- vgs[i-1,]
          ves[i,] <- ves[i-1,]
          ML1s[i,] <- ML1s[i-1,]
          ML0s[i,] <- ML0s[i-1,]
        }
        ps[i,] <- ps[i-1,]
      }
      else {
        if ( is.null(Z) ) {
          X <- cbind(X0[vids,,drop=FALSE],xs[i,vids])
          eig.R1 = emma.eigen.R.wo.Z(K[vids,vids],X)
        }
        else {
          vrows <- as.logical(rowSums(Z[,vids]))
          nr <- sum(vrows)
          X <- cbind(X0[vrows,,drop=FALSE],Z[vrows,vids]%*%t(xs[i,vids,drop=FALSE]))
          eig.R1 = emma.eigen.R.w.Z(Z[vrows,vids],K[vids,vids],X)
        }

        for(j in 1:g) {
          if ( nv == t ) {
            MLE <- emma.MLE(ys[j,],X,K,Z,ngrids,llim,ulim,esp,eig.L,eig.R1)
            if (!ponly) {
              ML1s[i,j] <- MLE$ML
              vgs[i,j] <- MLE$vg
              ves[i,j] <- MLE$ve
            }
            stats[i,j] <- 2*(MLE$ML-ML0[j])
          } else {
            if ( is.null(Z) ) {
              eig.L0 <- emma.eigen.L.wo.Z(K[vids,vids])
              MLE0 <- emma.MLE(ys[j,vids],X0[vids,,drop=FALSE],K[vids,vids],NULL,ngrids,llim,ulim,esp,eig.L0)
              MLE1 <- emma.MLE(ys[j,vids],X,K[vids,vids],NULL,ngrids,llim,ulim,esp,eig.L0)
            } else {
              if ( nr == n ) {
                MLE1 <- emma.MLE(ys[j,],X,K,Z,ngrids,llim,ulim,esp,eig.L)
              } else {
                eig.L0 <- emma.eigen.L.w.Z(Z[vrows,vids],K[vids,vids])
                MLE0 <- emma.MLE(ys[j,vrows],X0[vrows,,drop=FALSE],K[vids,vids],Z[vrows,vids],ngrids,llim,ulim,esp,eig.L0)
                MLE1 <- emma.MLE(ys[j,vrows],X,K[vids,vids],Z[vrows,vids],ngrids,llim,ulim,esp,eig.L0)
              }
            }
            if (!ponly) {
              ML1s[i,j] <- MLE1$ML
              ML0s[i,j] <- MLE0$ML
              vgs[i,j] <- MLE1$vg
              ves[i,j] <- MLE1$ve
            }
            stats[i,j] <- 2*(MLE1$ML-MLE0$ML)
          }
        }

        if ( ( nv == t ) && ( !ponly ) ) {
          ML0s[i,] <- ML0
        }
        ps[i,] <- pchisq(stats[i,],1,lower.tail=FALSE)
      }
    }
  }
  else {
    eig.L <- emma.eigen.L(Z,K)
    eig.R0 <- emma.eigen.R(Z,K,X0)

    for(i in 1:g) {
      vrows <- !is.na(ys[i,])
      if ( is.null(Z) ) {
        ML0[i] <- emma.MLE(ys[i,vrows],X0[vrows,,drop=FALSE],K[vrows,vrows],NULL,ngrids,llim,ulim,esp)$ML
      }
      else {
        vids <- colSums(Z[vrows,]>0)

        ML0[i] <- emma.MLE(ys[i,vrows],X0[vrows,,drop=FALSE],K[vids,vids],Z[vrows,vids],ngrids,llim,ulim,esp)$ML
      }
    }

    x.prev <- vector(length=0)

    # Initialize Progress Bar
    pb <- txtProgressBar(min = 0, initial = 0, max = m, style = 3)

    # Loops through each of the SNPs
    for(i in 1:m) {
      vids <- !is.na(xs[i,])
      nv <- sum(vids)
      xv <- xs[i,vids]

      # Update Progress Bar
      if (i %% 10 == 0) {
        Sys.sleep(0.001)
        setTxtProgressBar(pb, i)
      }

      if ( ( mean(xv) <= 0 ) || ( mean(xv) >= 1 ) ) {
        if (!ponly) {
          stats[i,] <- rep(NA,g)
          vgs[i,] <- rep(NA,g)
          ves[i,] <- rep(NA,g)
          ML1s[i,] <- rep(NA,g)
          ML0s[,i] <- rep(NA,g)
        }
        ps[i,] = rep(1,g)
      }
      else if ( identical(x.prev, xv) ) {
        if ( !ponly ) {
          stats[i,] <- stats[i-1,]
          vgs[i,] <- vgs[i-1,]
          ves[i,] <- ves[i-1,]
          ML1s[i,] <- ML1s[i-1,]
        }
        ps[i,] = ps[i-1,]
      }
      else {
        if ( is.null(Z) ) {
          X <- cbind(X0,xs[i,])
          if ( nv == t ) {
            eig.R1 = emma.eigen.R.wo.Z(K,X)
          }
        }
        else {
          vrows <- as.logical(rowSums(Z[,vids]))
          X <- cbind(X0,Z[,vids,drop=FALSE]%*%t(xs[i,vids,drop=FALSE]))
          if ( nv == t ) {
            eig.R1 = emma.eigen.R.w.Z(Z,K,X)
          }
        }

        for(j in 1:g) {
          vrows <- !is.na(ys[j,])
          if ( nv == t ) {
            nr <- sum(vrows)
            if ( is.null(Z) ) {
              if ( nr == n ) {
                MLE <- emma.MLE(ys[j,],X,K,NULL,ngrids,llim,ulim,esp,eig.L,eig.R1)
              }
              else {
                MLE <- emma.MLE(ys[j,vrows],X[vrows,],K[vrows,vrows],NULL,ngrids,llim,ulim,esp)
              }
            } else {
              if ( nr == n ) {
                MLE <- emma.MLE(ys[j,],X,K,Z,ngrids,llim,ulim,esp,eig.L,eig.R1)
              } else {
                vtids <- as.logical(colSums(Z[vrows,,drop=FALSE]))
                MLE <- emma.MLE(ys[j,vrows],X[vrows,],K[vtids,vtids],Z[vrows,vtids],ngrids,llim,ulim,esp)
              }
            }

            if (!ponly) {
              ML1s[i,j] <- MLE$ML
              vgs[i,j] <- MLE$vg
              ves[i,j] <- MLE$ve
            }
            stats[i,j] <- 2*(MLE$ML-ML0[j])
          }
          else {
            if ( is.null(Z) ) {
              vtids <- vrows & vids
              eig.L0 <- emma.eigen.L(NULL,K[vtids,vtids])
              MLE0 <- emma.MLE(ys[j,vtids],X0[vtids,,drop=FALSE],K[vtids,vtids],NULL,ngrids,llim,ulim,esp,eig.L0)
              MLE1 <- emma.MLE(ys[j,vtids],X[vtids,],K[vtids,vtids],NULL,ngrids,llim,ulim,esp,eig.L0)
            }
            else {
              vtids <- as.logical(colSums(Z[vrows,])) & vids
              vtrows <- vrows & as.logical(rowSums(Z[,vids]))
              eig.L0 <- emma.eigen.L(Z[vtrows,vtids],K[vtids,vtids])
              MLE0 <- emma.MLE(ys[j,vtrows],X0[vtrows,,drop=FALSE],K[vtids,vtids],Z[vtrows,vtids],ngrids,llim,ulim,esp,eig.L0)
              MLE1 <- emma.MLE(ys[j,vtrows],X[vtrows,],K[vtids,vtids],Z[vtrows,vtids],ngrids,llim,ulim,esp,eig.L0)
            }
            if (!ponly) {
              ML1s[i,j] <- MLE1$ML
              vgs[i,j] <- MLE1$vg
              ves[i,j] <- MLE1$ve
              ML0s[i,j] <- MLE0$ML
            }
            stats[i,j] <- 2*(MLE1$ML-MLE0$ML)
          }
        }
        if ( ( nv == t ) && ( !ponly ) ) {
          ML0s[i,] <- ML0
        }
        ps[i,] <- pchisq(stats[i,],1,lower.tail=FALSE)
      }
    }
  }
  if ( ponly ) {
    return (ps)
  }
  else {
    return (list(ps=ps,ML1s=ML1s,ML0s=ML0s,stats=stats,vgs=vgs,ves=ves))
  }
}
