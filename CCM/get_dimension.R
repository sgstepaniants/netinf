library(rEDM)

# Number of Nodes
n = 5;
# Time Vector
time = 1:500
# Boundary Conditions
bc = "circ"

Es <- 1:(2*n+1)
rho_ave <- numeric(2*n+1)
trials <- 100
for (i in 1:trials) {
  # Generate Nearest-Neighbor Spring Mass Model Data
  ts <- nncoupled_data(n, time, smallrandpfn, zerovfn, constmfn, constkfn, bc=bc)
  
  # portion of the data used to create reconstruction
  lib <- c(1, 100)
  # the reconstruction will be used to make forecasts on this portion
  pred <- c(201, 500)
  
  # get a data frame for the model parameters and forecast statistics
  simplex_output <- simplex(ts, lib, pred, E=Es)
  rho_ave <- rho_ave + simplex_output$rho
}
rho_ave <- rho_ave / trials

par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0))  # set up margins for plotting
plot(Es, rho_ave, type = "l", xlab = "Embedding Dimension (E)", 
     ylab = "Forecast Skill (rho)", main=paste("n = ", + n))



