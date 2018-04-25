library(rEDM)
library(lattice)

# Number of Nodes
n <- 4
E <- 4
# Time Vector
tmax <- 100
time <- 1:tmax
# Boundary Conditions
bc <- "circ"

trials <- 100
libs <- seq(10, tmax-10, by = 10)
noise <- 0.05

# position function that places nodes 2 and 4 at equilibrium
pinpfn <- function(n, bc) {
  return(c(runif(1, -0.25, 0.25), 0.25, runif(1, 0.25, 0.75), 0.75))
}
# mass function that makes masses of nodes 2 and 4 very large
pinmfn <- function(n) {
  return(c(1, 10000, 1, 10000))
}

# Generate Nearest-Neighbor Spring Mass Model Adjacency Matrix
names <- paste("pos", 1:n, sep="")
data_func <- function() {
  nncoupled_model <- nncoupled_data(n, time, pinpfn, zerovfn,
                                      pinmfn, constkfn, bc=bc)
  nncoupled_model <- nncoupled_model + cbind(0,matrix(rnorm(n*tmax,mean=0,sd=noise),tmax,n))
  return(nncoupled_model)
}

adj <- get_ccm_means(data_func, E, names, lib_sizes=libs, trials=trials)
adj_comp <- matrix(0, length(names), length(names))
adj_opt_rho <- adj_opt_dim
max_rho <- apply(adj, c(1, 2), max)
# Search through rho values and find the one that represents 95%
# of the correlation maximum for each graph. Get the embedding dimension
# of this point also.
for (i in 1:length(names)) {
  for (j in 1:length(names)) {
    idx <- which(adj[i, j,] > 0.95 * max_rho[i, j])[1]
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
          at=seq(-1, 1, length.out=120),
          panel=myPanel)

