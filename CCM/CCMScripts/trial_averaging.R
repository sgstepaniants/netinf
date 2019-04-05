library(rEDM)
library(lattice)
#source('../generate_data.R')
#source('../ccm_helper.R')

# Number of Nodes
n <- 4
E <- 4
# Time Vector
tmax <- 200
time <- 1:tmax
# Boundary Conditions
bc <- "circ"

trials <- 100
libs <- seq(10, tmax-10, by = 10)
noise <- 0

# Generate Nearest-Neighbor Spring Mass Model Adjacency Matrix
names <- paste("pos", 1:n, sep="")
data_func <- function() {
  nncoupled_model <- nncoupled_data(n, time, randpfn, zerovfn,
                                      largerandmfn, largerandkfn, bc=bc)
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
