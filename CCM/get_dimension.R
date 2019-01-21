library(rEDM)

# Number of Nodes
n <- 30;
# Time Vector
tmax <- 30
deltat = tmax / 200;
time <- seq(1, tmax, by=deltat)
# Boundary Conditions
bc <- "fixed"

Es <- 1:(2*n+1)
rho_ave <- numeric(2*n+1)

libs <- seq(10, tmax / deltat - 10, by=20)
noise <- 0.01

p <- 0.5
k_const <- 0.01

pert <- 1
onepfn <- function(n, bc) {
  p <- unifpfn(n, bc);
  i <- sample(1 : n, 1)
  p[i] <- p[i] + pert * rnorm(1)
  return(p)
}

damping <- 0.25
constcfn <- function(n) {
  return(rep(damping, n))
}

data_func <- function() {
  # Construct the adjacency matrix
  A <- symmetric_erdos_reyni(n, p)
  
  # Construct the spring constants matrix
  K <- tridiag(n + 2, corners=FALSE)
  K[2 : (n + 1), 2 : (n + 1)] <- A
  K <- k_const * K
  
  nncoupled_model <- nncoupled_data(n, time, K, onepfn, zerovfn,
                                    constmfn, constcfn, bc=bc, pert=0);
  nncoupled_model <- nncoupled_model + cbind(0,matrix(rnorm(n*length(time),mean=0,sd=noise),length(time),n))
  
  matplot(nncoupled_model[, !(colnames(nncoupled_model) %in% "time")], col=rainbow(n), type="l", lty=1)
  legend("topright", legend=names, col=rainbow(n), lty=1, cex=0.8)
  return(nncoupled_model)
}

trials <- 10
for (i in 1:trials) {
  print(i)
  # Generate Nearest-Neighbor Spring Mass Model Data
  ts <- data_func()
  
  # portion of the data used to create reconstruction
  div <- as.integer(length(time) / 3)
  lib <- c(1, div)
  # the reconstruction will be used to make forecasts on this portion
  pred <- c(div + 1, length(time))
  
  # get a data frame for the model parameters and forecast statistics
  simplex_output <- simplex(ts, lib, pred, E=Es)
  rho_ave <- rho_ave + simplex_output$rho
}
rho_ave <- rho_ave / trials

par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0))  # set up margins for plotting
plot(Es, rho_ave, type = "l", xlab = "Embedding Dimension (E)", 
     ylab = "Forecast Skill (rho)", main=paste("n = ", n, ", p = ", p))
