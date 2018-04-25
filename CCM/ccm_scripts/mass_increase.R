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

adj <- get_ccm_means(data_func, E, names, lib_sizes=libs, trials=trials)
adj_comp <- matrix(0, length(names), length(names))
adj_opt_rho <- adj_comp
max_rho <- apply(adj, c(1, 2), max)
# Search through rho values and find the one that represents 95%
# of the correlation maximum for each graph. Get the embedding dimension
# of this point also.
for (i in 1:length(names)) {
  for (j in 1:length(names)) {
    idx <- which(adj[i, j,] >= 0.95 * max_rho[i, j])[1]
    adj_comp[i, j] <- libs[idx]
    adj_opt_rho[i, j] <- adj[i, j, idx]
  }
}

adj_direc <- adj_opt_rho - t(adj_opt_rho)

# complexity is optimal embedding dimension
rgb.palette <- colorRampPalette(c("blue", "red"), space = "rgb")
myPanel <- function(x, y, z, ...) {
  panel.levelplot(x,y,z,...)
  panel.text(x, y, round(z,3))
}

levelplot(t(adj_comp), main="NNCoupled Complexity",
          xlab="", ylab="", ylim=c(length(names) + 0.5, 0.5),
          col.regions=rgb.palette(120),
          at=seq(0, max(libs), length.out=120),
          panel=myPanel)
levelplot(t(adj_direc), main="NNCoupled Directionality",
          xlab="", ylab="", ylim=c(length(names) + 0.5, 0.5),
          col.regions=rgb.palette(120),
          at=seq(-1, 1, length.out=120),
          panel=myPanel)
levelplot(t(adj_opt_rho), main="NNCoupled Correlation",
          xlab="", ylab="", ylim=c(length(names) + 0.5, 0.5),
          col.regions=rgb.palette(120),
          at=seq(0, 1, length.out=120),
          panel=myPanel)
