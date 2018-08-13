emma.kinship <- function(snps, method="additive", use="all") {

  n0 <- sum(snps==0,na.rm=TRUE)
  nh <- sum(snps==0.5,na.rm=TRUE)
  n1 <- sum(snps==1,na.rm=TRUE)
  nNA <- sum(is.na(snps))

  stopifnot(n0+nh+n1+nNA == length(snps))

  if ( method == "dominant" ) {
    flags <- matrix(as.double(rowMeans(snps,na.rm=TRUE) > 0.5),nrow(snps),ncol(snps))
    snps[!is.na(snps) & (snps == 0.5)] <- flags[!is.na(snps) & (snps == 0.5)]
  }
  else if ( method == "recessive" ) {
    flags <- matrix(as.double(rowMeans(snps,na.rm=TRUE) < 0.5),nrow(snps),ncol(snps))
    snps[!is.na(snps) & (snps == 0.5)] <- flags[!is.na(snps) & (snps == 0.5)]
  }
  else if ( ( method == "additive" ) && ( nh > 0 ) ) {
    dsnps <- snps
    rsnps <- snps
    flags <- matrix(as.double(rowMeans(snps,na.rm=TRUE) > 0.5),nrow(snps),ncol(snps))
    dsnps[!is.na(snps) & (snps==0.5)] <- flags[!is.na(snps) & (snps==0.5)]
    flags <- matrix(as.double(rowMeans(snps,na.rm=TRUE) < 0.5),nrow(snps),ncol(snps))
    rsnps[!is.na(snps) & (snps==0.5)] <- flags[!is.na(snps) & (snps==0.5)]
    snps <- rbind(dsnps,rsnps)
  }

  if ( use == "all" ) {
    mafs <- matrix(rowMeans(snps,na.rm=TRUE),nrow(snps),ncol(snps))
    snps[is.na(snps)] <- mafs[is.na(snps)]
  }
  else if ( use == "complete.obs" ) {
    snps <- snps[rowSums(is.na(snps))==0,]
  }

  n <- ncol(snps)
  K <- matrix(nrow=n,ncol=n)
  diag(K) <- 1

  # Initialize Progress Bar
  pb <- txtProgressBar(min = 0, initial = 0, max = n, style = 3)

  for(i in 2:n) {
    for(j in 1:(i-1)) {
      x <- snps[,i]*snps[,j] + (1-snps[,i])*(1-snps[,j])
      K[i,j] <- sum(x,na.rm=TRUE)/sum(!is.na(x))
      K[j,i] <- K[i,j]
    }

    # Update Progress Bar
    Sys.sleep(0.001)
    setTxtProgressBar(pb, i)
  }
  return(K)
}
