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
