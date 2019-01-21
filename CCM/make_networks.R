tridiag <- function(n, corners=TRUE) {
  out <- matrix(0, n, n)
  indx <- seq.int(n - 1)
  out[cbind(indx+1,indx)] <- rep(1, n - 1)
  out[cbind(indx,indx+1)] <- rep(1, n - 1)
  
  if (corners) {
    out[1, n] <- 1
    out[n, 1] <- 1
  }
  
  return(out)
}

erdos_reyni <- function(n, p) {
  out <- matrix(0, n, n)
  mat <- matrix(runif(n * n), n, n)

  # if >= 1-p, add connection
  out[mat >= 1 - p] <- 1
  return(out)
}

symmetric_erdos_reyni <- function(n, p, zero_diag=TRUE) {
  er <- erdos_reyni(n, p)
  er[lower.tri(er, diag=zero_diag)] <- 0
  
  out <- matrix(0, n, n)
  out <- er + t(er)
  return(out)
}
