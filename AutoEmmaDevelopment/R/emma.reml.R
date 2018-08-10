# EMMA using REML
emma.REML.t <- function(ys, xs, K, Z=NULL, X0 = NULL, ngrids=100, llim=-10, ulim=10, esp=1e-10, ponly = FALSE) {
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

  stopifnot(nrow(K) == t)
  stopifnot(ncol(K) == t)
  stopifnot(nrow(X0) == n)

  if ( !ponly ) {
    REMLs <- matrix(nrow=m,ncol=g)
    vgs <- matrix(nrow=m,ncol=g)
    ves <- matrix(nrow=m,ncol=g)
  }
  dfs <- matrix(nrow=m,ncol=g)
  stats <- matrix(nrow=m,ncol=g)
  ps <- matrix(nrow=m,ncol=g)

  if ( sum(is.na(ys)) == 0 ) {
    eig.L <- emma.eigen.L(Z,K)

    x.prev <- vector(length=0)

    # Initialize Progress Bar
    pb <- txtProgressBar(min = 0, initial = 0, max = m, style = 3)

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
        if ( !ponly ) {
          vgs[i,] <- rep(NA,g)
          ves[i,] <- rep(NA,g)
          dfs[i,] <- rep(NA,g)
          REMLs[i,] <- rep(NA,g)
          stats[i,] <- rep(NA,g)
        }
        ps[i,] = rep(1,g)

      }
      else if ( identical(x.prev, xv) ) {
        if ( !ponly ) {
          vgs[i,] <- vgs[i-1,]
          ves[i,] <- ves[i-1,]
          dfs[i,] <- dfs[i-1,]
          REMLs[i,] <- REMLs[i-1,]
          stats[i,] <- stats[i-1,]
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
          X <- cbind(X0[vrows,,drop=FALSE],Z[vrows,vids,drop=FALSE]%*%t(xs[i,vids,drop=FALSE]))
          eig.R1 = emma.eigen.R.w.Z(Z[vrows,vids],K[vids,vids],X)
        }

        for(j in 1:g) {
          if ( nv == t ) {
            REMLE <- emma.REMLE(ys[j,],X,K,Z,ngrids,llim,ulim,esp,eig.R1)
            if ( is.null(Z) ) {
              U <- eig.L$vectors * matrix(sqrt(1/(eig.L$values+REMLE$delta)),t,t,byrow=TRUE)
              dfs[i,j] <- nv - q1
            }
            else {
              U <- eig.L$vectors * matrix(c(sqrt(1/(eig.L$values+REMLE$delta)),rep(sqrt(1/REMLE$delta),n-t)),n,n,byrow=TRUE)
              dfs[i,j] <- n - q1
            }
            yt <- crossprod(U,ys[j,])
            Xt <- crossprod(U,X)
            iXX <- solve(crossprod(Xt,Xt))
            beta <- iXX%*%crossprod(Xt,yt)

            if ( !ponly ) {
              vgs[i,j] <- REMLE$vg
              ves[i,j] <- REMLE$ve
              REMLs[i,j] <- REMLE$REML
            }
            stats[i,j] <- beta[q1]/sqrt(iXX[q1,q1]*REMLE$vg)
          }
          else {
            if ( is.null(Z) ) {
              eig.L0 <- emma.eigen.L.wo.Z(K[vids,vids])
              nr <- sum(vids)
              yv <- ys[j,vids]
              REMLE <- emma.REMLE(yv,X,K[vids,vids,drop=FALSE],NULL,ngrids,llim,ulim,esp,eig.R1)
              U <- eig.L0$vectors * matrix(sqrt(1/(eig.L0$values+REMLE$delta)),nr,nr,byrow=TRUE)
              dfs[i,j] <- nr - q1
            }
            else {
              eig.L0 <- emma.eigen.L.w.Z(Z[vrows,vids,drop=FALSE],K[vids,vids])
              yv <- ys[j,vrows]
              nr <- sum(vrows)
              tv <- sum(vids)
              REMLE <- emma.REMLE(yv,X,K[vids,vids,drop=FALSE],Z[vrows,vids,drop=FALSE],ngrids,llim,ulim,esp,eig.R1)
              U <- eig.L0$vectors * matrix(c(sqrt(1/(eig.L0$values+REMLE$delta)),rep(sqrt(1/REMLE$delta),nr-tv)),nr,nr,byrow=TRUE)
              dfs[i,j] <- nr - q1
            }
            yt <- crossprod(U,yv)
            Xt <- crossprod(U,X)
            iXX <- solve(crossprod(Xt,Xt))
            beta <- iXX%*%crossprod(Xt,yt)
            if (!ponly) {
              vgs[i,j] <- REMLE$vg
              ves[i,j] <- REMLE$ve
              REMLs[i,j] <- REMLE$REML
            }
            stats[i,j] <- beta[q1]/sqrt(iXX[q1,q1]*REMLE$vg)
          }
        }
        ps[i,] <- 2*pt(abs(stats[i,]),dfs[i,],lower.tail=FALSE)
      }
    }
  }
  else {
    eig.L <- emma.eigen.L(Z,K)
    eig.R0 <- emma.eigen.R(Z,K,X0)

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
          vgs[i,] <- rep(NA,g)
          ves[i,] <- rep(NA,g)
          REMLs[i,] <- rep(NA,g)
          dfs[i,] <- rep(NA,g)
        }
        ps[i,] = rep(1,g)
      }
      else if ( identical(x.prev, xv) ) {
        if ( !ponly ) {
          stats[i,] <- stats[i-1,]
          vgs[i,] <- vgs[i-1,]
          ves[i,] <- ves[i-1,]
          REMLs[i,] <- REMLs[i-1,]
          dfs[i,] <- dfs[i-1,]
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
          vrows <- as.logical(rowSums(Z[,vids,drop=FALSE]))
          X <- cbind(X0,Z[,vids,drop=FALSE]%*%t(xs[i,vids,drop=FALSE]))
          if ( nv == t ) {
            eig.R1 = emma.eigen.R.w.Z(Z,K,X)
          }
        }

        for(j in 1:g) {
          vrows <- !is.na(ys[j,])
          if ( nv == t ) {
            yv <- ys[j,vrows]
            nr <- sum(vrows)
            if ( is.null(Z) ) {
              if ( nr == n ) {
                REMLE <- emma.REMLE(yv,X,K,NULL,ngrids,llim,ulim,esp,eig.R1)
                U <- eig.L$vectors * matrix(sqrt(1/(eig.L$values+REMLE$delta)),n,n,byrow=TRUE)
              }
              else {
                eig.L0 <- emma.eigen.L.wo.Z(K[vrows,vrows,drop=FALSE])
                REMLE <- emma.REMLE(yv,X[vrows,,drop=FALSE],K[vrows,vrows,drop=FALSE],NULL,ngrids,llim,ulim,esp)
                U <- eig.L0$vectors * matrix(sqrt(1/(eig.L0$values+REMLE$delta)),nr,nr,byrow=TRUE)
              }
              dfs[i,j] <- nr-q1
            }
            else {
              if ( nr == n ) {
                REMLE <- emma.REMLE(yv,X,K,Z,ngrids,llim,ulim,esp,eig.R1)
                U <- eig.L$vectors * matrix(c(sqrt(1/(eig.L$values+REMLE$delta)),rep(sqrt(1/REMLE$delta),n-t)),n,n,byrow=TRUE)
              }
              else {
                vtids <- as.logical(colSums(Z[vrows,,drop=FALSE]))
                eig.L0 <- emma.eigen.L.w.Z(Z[vrows,vtids,drop=FALSE],K[vtids,vtids,drop=FALSE])
                REMLE <- emma.REMLE(yv,X[vrows,,drop=FALSE],K[vtids,vtids,drop=FALSE],Z[vrows,vtids,drop=FALSE],ngrids,llim,ulim,esp)
                U <- eig.L0$vectors * matrix(c(sqrt(1/(eig.L0$values+REMLE$delta)),rep(sqrt(1/REMLE$delta),nr-sum(vtids))),nr,nr,byrow=TRUE)
              }
              dfs[i,j] <- nr-q1
            }

            yt <- crossprod(U,yv)
            Xt <- crossprod(U,X[vrows,,drop=FALSE])
            iXX <- solve(crossprod(Xt,Xt))
            beta <- iXX%*%crossprod(Xt,yt)
            if ( !ponly ) {
              vgs[i,j] <- REMLE$vg
              ves[i,j] <- REMLE$ve
              REMLs[i,j] <- REMLE$REML
            }
            stats[i,j] <- beta[q1]/sqrt(iXX[q1,q1]*REMLE$vg)
          }
          else {
            if ( is.null(Z) ) {
              vtids <- vrows & vids
              eig.L0 <- emma.eigen.L.wo.Z(K[vtids,vtids,drop=FALSE])
              yv <- ys[j,vtids]
              nr <- sum(vtids)
              REMLE <- emma.REMLE(yv,X[vtids,,drop=FALSE],K[vtids,vtids,drop=FALSE],NULL,ngrids,llim,ulim,esp)
              U <- eig.L0$vectors * matrix(sqrt(1/(eig.L0$values+REMLE$delta)),nr,nr,byrow=TRUE)
              Xt <- crossprod(U,X[vtids,,drop=FALSE])
              dfs[i,j] <- nr-q1
            }
            else {
              vtids <- as.logical(colSums(Z[vrows,,drop=FALSE])) & vids
              vtrows <- vrows & as.logical(rowSums(Z[,vids,drop=FALSE]))
              eig.L0 <- emma.eigen.L.w.Z(Z[vtrows,vtids,drop=FALSE],K[vtids,vtids,drop=FALSE])
              yv <- ys[j,vtrows]
              nr <- sum(vtrows)
              REMLE <- emma.REMLE(yv,X[vtrows,,drop=FALSE],K[vtids,vtids,drop=FALSE],Z[vtrows,vtids,drop=FALSE],ngrids,llim,ulim,esp)
              U <- eig.L0$vectors * matrix(c(sqrt(1/(eig.L0$values+REMLE$delta)),rep(sqrt(1/REMLE$delta),nr-sum(vtids))),nr,nr,byrow=TRUE)
              Xt <- crossprod(U,X[vtrows,,drop=FALSE])
              dfs[i,j] <- nr-q1
            }
            yt <- crossprod(U,yv)
            iXX <- solve(crossprod(Xt,Xt))
            beta <- iXX%*%crossprod(Xt,yt)
            if ( !ponly ) {
              vgs[i,j] <- REMLE$vg
              ves[i,j] <- REMLE$ve
              REMLs[i,j] <- REMLE$REML
            }
            stats[i,j] <- beta[q1]/sqrt(iXX[q1,q1]*REMLE$vg)

          }
        }
        ps[i,] <- 2*pt(abs(stats[i,]),dfs[i,],lower.tail=FALSE)
      }
    }
  }
  if ( ponly ) {
    return (ps)
  }
  else {
    return (list(ps=ps,REMLs=REMLs,stats=stats,dfs=dfs,vgs=vgs,ves=ves))
  }
}
