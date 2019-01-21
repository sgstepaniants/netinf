library(rEDM)
library(lattice)
#source('../generate_data.R')
#source('../ccm_helper.R')

# Number of Nodes
n <- 3
E <- 3

# Information to construct Hankel matrix
stacks <- 500 # times to stack data in Hankel matrix
rank <- n # rank to truncate when taking SVD of Hankel matrix

# Time Vector
deltat = 0.005;
tmax <- 3
time <- seq(1, tmax + deltat * (stacks - 1), by=deltat)
# Boundary Conditions
bc <- "fixed"

trials <- 20
libs <- seq(10, tmax / deltat - 10, by=50)
noise <- 0

# Generate Nearest-Neighbor Spring Mass Model Adjacency Matrix
names <- paste("u", 1 : rank, sep="")
data_func <- function() {
  nncoupled_model <- nncoupled_data(n, time, randpfn, zerovfn,
                                      constmfn, constkfn, zerocfn, bc=bc, pert=0)
  #nncoupled_model <- nncoupled_model + cbind(0,matrix(rnorm(n*tmax,mean=0,sd=noise),tmax,n))
  
  matplot(nncoupled_model[, !(colnames(nncoupled_model) %in% "time")], type="l")
  data_hankel <- hankel(nncoupled_model, stacks);
  return(data_hankel);
}

ccm_rho_hankel <- get_ccm_rho_hankel(data_func, E, n, names, lib_sizes=libs, trials=trials);
ccm_rho_graphs <- ccm_rho_hankel$ccm_rho_graphs;
svd_data <- ccm_rho_hankel$svd_data;
adj <- abs(get_adj_hankel(ccm_rho_graphs, svd_data));
adj_pool <- normalize(pool(adj, stacks, stacks, "sum"));

rgb.palette <- colorRampPalette(c("blue", "red"), space = "rgb")
myPanel <- function(x, y, z, ...) {
  panel.levelplot(x,y,z,...)
  panel.text(x, y, round(z,3))
}

levelplot(adj_pool, main="NNCoupled Correlation",
          xlab="", ylab="", ylim=c(length(names) + 0.5, 0.5),
          col.regions=rgb.palette(120),
          at=seq(0, 1, length.out=120),
          panel=myPanel)
