library(rEDM)
library(lattice)
#source('../generate_data.R')
#source('../ccm_helper.R')

# Number of Nodes
n <- 4
E <- 4
# Time Vector
tmax <- 20
time <- seq(1, tmax, by=0.1)
# Boundary Conditions
bc <- "circ"

trials <- 10
libs <- seq(10, length(time)-10, by=20)
noise <- 0.01

randpfn <- function(n, bc) {
  return(c(0.1, 0.25, 0.5, 0.75))
}

constmfn <- function(n) {
  return(c(100, 100, 100, 100))
}

constkfn <- function(n, bc) {
  k <- rep(1, n)  # circ
  if (bc == "fixed") {
    k <- rep(50, n+1)
  } else if (bc == "free") {
    k <- rep(50, n-1)
  }
  return(k)
}

# Generate Nearest-Neighbor Spring Mass Model Adjacency Matrix
names <- paste("pos", 1:n, sep="")
data_func <- function() {
  nncoupled_model <- nncoupled_data(n, time, randpfn, zerovfn,
                                      constmfn, constkfn, zerocfn, bc=bc, pert=0)
  nncoupled_model <- nncoupled_model + cbind(0,matrix(rnorm(n*length(time),mean=0,sd=noise),length(time),n))
  matplot(nncoupled_model[, !(colnames(nncoupled_model) %in% "time")], type="l")
  return(nncoupled_model)
}

ccm_rho_graphs <- get_ccm_rho(data_func, E, names, lib_sizes=libs, trials=trials);
adj <- get_adj(ccm_rho_graphs);


rgb.palette <- colorRampPalette(c("blue", "red"), space = "rgb")
myPanel <- function(x, y, z, ...) {
  panel.levelplot(x,y,z,...)
  panel.text(x, y, round(z,3))
}

levelplot(adj, main="NNCoupled Correlation",
          xlab="", ylab="", ylim=c(length(names) + 0.5, 0.5),
          col.regions=rgb.palette(120),
          at=seq(0, 1, length.out=120),
          panel=myPanel)
