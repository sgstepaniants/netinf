library(rEDM)
library(lattice)

# Number of Nodes
n <- 10
E <- 10
# Time Vector
tmax <- 10
deltat = tmax / 200;
time <- seq(1, tmax, by=deltat)
# Boundary Conditions
bc <- "fixed"

trials <- 100
libs <- seq(10, tmax / deltat - 10, by=20)
noise <- 0.01

onepfn <- function(n, bc) {
  p <- unifpfn(n, bc);
  i <- sample(1 : n, 1)
  p[i] <- p[i] - 0.5
  return(p)
}

k <- 0.1
constkfn <- function(n, bc) {
  return(rep(k, n + 1))
}

damping <- 0.5
constcfn <- function(n) {
  return(rep(damping, n))
}

# Generate Nearest-Neighbor Spring Mass Model Adjacency Matrix
names <- paste("pos", 1:n, sep="")
data_func <- function() {
  nncoupled_model <- nncoupled_data(n, time, onepfn, zerovfn,
                                    constmfn, constkfn, constcfn, bc=bc, pert=0);
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
