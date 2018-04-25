library(rEDM)

# Number of Nodes
n <- 4
E <- 41
# Time Vector
tmax <- 100
time <- 1:tmax
# Boundary Conditions
bc <- "circ"

trials <- 10
libs <- seq(10, tmax-10, by = 10)

# position function that places nodes 2 and 4 at equilibrium
pinpfn <- function(n, bc) {
  return(c(runif(1, -0.25, 0.25), 0.25, runif(1, 0.25, 0.75), 0.75))
}
# mass function that makes masses of nodes 2 and 4 very large
pinmfn <- function(n) {
  return(c(1, 10000, 1, 10000))
}

# Generate Nearest-Neighbor Spring Mass Model Data
a_xmap_b_rho <- numeric(length(libs))
b_xmap_a_rho <- a_xmap_b_rho
for (i in 1:trials) {
  nncoupled_model <- nncoupled_data(n, time, smallrandpfn, zerovfn,
                                    constmfn, constkfn, bc=bc)
  nncoupled_model <- nncoupled_model + cbind(0,matrix(rnorm(n*tmax,mean=0,sd=0.1),tmax,n))
  
  a_xmap_b <- ccm(nncoupled_model, E = E, lib_column = "pos1", 
                  target_column = "pos3", lib_sizes = libs,
                  num_samples = 100, random_libs = TRUE, replace = TRUE);
  b_xmap_a <- ccm(nncoupled_model, E = E, lib_column = "pos3",
                  target_column = "pos1", lib_sizes = libs,
                  num_samples = 100, random_libs = TRUE, replace = TRUE);
  
  a_xmap_b_means <- ccm_means(a_xmap_b);
  b_xmap_a_means <- ccm_means(b_xmap_a);
  a_xmap_b_rho <- a_xmap_b_rho + a_xmap_b_means$rho
  b_xmap_a_rho <- b_xmap_a_rho + b_xmap_a_means$rho
  print(i)
}
a_xmap_b_rho <- a_xmap_b_rho / trials
b_xmap_a_rho <- b_xmap_a_rho / trials

par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0))
ab <- pmax(0, a_xmap_b_rho);
ba <- pmax(0, b_xmap_a_rho);

plot(a_xmap_b_means$lib_size, ab, type = "l", col = "red", xlab = "Library Size", 
     ylab = "Cross Map Skill (rho)", ylim = c(0, 1), main = "a = 1, b = 3")
lines(b_xmap_a_means$lib_size, ba, col = "blue")
legend(x = "topleft", legend = c("a xmap b", "b xmap a"),
       col = c("red", "blue"), lwd = 1, bty = "n", inset = 0.02, cex = 0.8)
