# Initialization functions for spring mass nearest neighbor model
randpfn <- function(n, bc) {
  p <- runif(n, -0.5, 0.5)
  return(p)
}

smallrandpfn <- function(n, bc) {
  p <- seq(from=-0.5, to=(n-1)/n - 0.5, by=1/n) + rnorm(n, mean=0, sd=1/(2*n)) # circ
  if (bc == "fixed") {
    p <- seq(from=1/(n+1) - 0.5, to=n/(n+1) - 0.5, by=1/(n+1)) + rnorm(n, mean=0, sd=1/(2*(n+1)))
  } else if (bc == "free") {
    p <- seq(from=-0.5, to=0.5, by=1/(n-1)) + rnorm(n, mean=0, sd=1/(2*(n-1)))
  }
  return(p)
}

singlepfn <- function(n, bc) {
  p <- c(-1/n, seq(from=1/n, to=(n-1)/n, by=1/n))  # circ
  if (bc == "fixed") {
    p <- c(0, seq(from=2/(n+1), to=n/(n+1), by=1/(n+1)))
  } else if (bc == "free") {
    p <- c(-1/(n-1), seq(from=1/(n-1), to=1, by=1/(n-1)))
  }
  return(p)
}

unifpfn <- function(n, bc) {
  p <- seq(from=0, to=(n-1), by=1) / n - (n-1)/(2*n)  # circ
  if (bc == "fixed") {
    p <- seq(from=1, to=n, by=1) / (n+1) - 0.5
  } else if (bc == "free") {
    if (n == 1) {
      return(0)
    }
    p <- seq(from=0, to=n-1, by=1) / (n-1) - 0.5
  }
  return(p)
}

zerovfn <- function(n) {
  return(rep(0, n))
}

constmfn <- function(n) {
  return(rep(1, n))
}

randmfn <- function(n) {
  return(runif(n, 0, 1))
}

largerandmfn <- function(n) {
  return(runif(n, 0, 100))
}

constkfn <- function(n, bc) {
  k <- rep(1, n)  # circ
  if (bc == "fixed") {
    k <- rep(1, n+1)
  } else if (bc == "free") {
    k <- rep(1, n-1)
  }
  return(k)
}

randkfn <- function(n, bc) {
  k <- runif(n, 0, 1) # circ
  if (bc == "fixed") {
    k <- runif(n+1, 0, 1)
  } else if (bc == "free") {
    k <- runif(n-1, 0, 1)
  }
  return(k)
}

largerandkfn <- function(n, bc) {
  k <- runif(n, 0, 100) # circ
  if (bc == "fixed") {
    k <- runif(n+1, 0, 100)
  } else if (bc == "free") {
    k <- runif(n-1, 0, 100)
  }
  return(k)
}

constcfn <- function(n) {
  return(rep(0.3, n))
}

zerocfn <- function(n) {
  return(rep(0, n))
}






# Initialization function for Kuramoto model
randicfn <- function(n) {
  ic <- 2*pi*runif(n, 0, 1);  # uniform [0, 2pi]
  return(ic)
}

randwfn <- function(n) {
  w <- 2*runif(n, 0, 1)  # uniform [0, 2]
  return(w)
}
