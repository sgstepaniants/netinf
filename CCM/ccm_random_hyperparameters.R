library(rEDM)
library(lattice)

# Experiment Name
exp_name <- "FirstExp"

# Information For Each Experiment
# Number of Nodes
n <- 5
E <- 5
# Time Vector
tmax <- 10
deltat = tmax / 200;
time <- seq(0, tmax, by=deltat)
# Boundary Conditions
bc <- "fixed"

# Number of trials in each experiment
trials <- 10
libs <- seq(10, tmax / deltat - 10, by=20)
noise <- 0.01

# Try out different probabilities
probs <- seq(0, 1, by=0.1)

# Construct the adjacency matrix
numMats <- 3
mats <- array(0, dim = c(length(probs), numMats, n, n))
for (i in 1 : length(probs)) {
  p <- probs[i]
  for (j in 1 : numMats) {
    mats[i, j,,] <- symmetric_erdos_reyni(n, p)
  }
}

# Node position initialization function
pert <- 1
onepfn <- function(n, bc) {
  pos <- unifpfn(n, bc);
  i <- sample(1 : n, 1)
  pos[i] <- pos[i] + pert * rnorm(1)
  return(pos)
}

# Try out different spring constants
k_consts <- seq(0.01, 0.1, by=0.01)

# Try out different damping coefficients
damping_coeffs <- seq(0.25, 0.5, by=0.25)

# Arrays to save the adjacency matrices and RMSE of each experiment
adj_results <- array(0, dim = c(length(probs), length(k_consts),
                                length(damping_coeffs),
                                numMats, n, n))
rmse_results <- array(0, dim = c(length(probs), length(k_consts),
                                 length(damping_coeffs), numMats))

out_file <- paste("~/netinf/CCM/experiments/", exp_name, ".rda", sep="")

# Run CCM experiments with different probabilities of connection,
# random connectivity matrices, spring constants,
# and damping coefficients.
for (l in 1 : length(probs)) {
  for (i in 1 : length(k_consts)) {
    for (j in 1 : length(damping_coeffs)) {
      for (m in 1 : numMats) {
        k_const <- k_consts[i]
        damping_coeff <- damping_coeffs[j]
        
        sprintf("p_num: %i, k_num: %i, damp_num: %i, mat_num: %i", l, i, j, m)
        truth <- mats[l, m,,]
        
        
        # Make the matrix of spring coefficients
        K <- tridiag(n + 2, corners=FALSE)
        K[2 : (n + 1), 2 : (n + 1)] <- truth
        K <- k_const * K
        
        # Make the damping coefficients function
        constcfn <- function(n) {
          return(rep(damping_coeff, n))
        }
        
        # Generate Nearest-Neighbor Spring Mass Model Adjacency Matrix
        names <- paste("pos", 1:n, sep="")
        data_func <- function() {
          nncoupled_model <- nncoupled_data(n, time, K, onepfn, zerovfn,
                                            constmfn, constcfn, bc=bc, pert=0);
          nncoupled_model <- nncoupled_model + cbind(0,matrix(rnorm(n*length(time),mean=0,sd=noise),length(time),n))
          
          matplot(nncoupled_model[, !(colnames(nncoupled_model) %in% "time")], col=rainbow(n), type="l", lty=1)
          legend("topright", legend=names, col=rainbow(n), lty=1, cex=0.8)
          return(nncoupled_model)
        }
        
        ccm_rho_graphs <- get_ccm_rho(data_func, E, names, lib_sizes=libs, trials=trials);
        adj_mats <- get_adj(ccm_rho_graphs);
        adj <- apply(adj_mats, c(1, 2), mean)
        
        # Zero out the diagonal of the resulting
        # adjacency matrix.
        diag(adj) <- 0
        
        # Save the results of each experiment
        adj_results[l, i, j, m,,] <- adj
        rmse_results[l, i, j, m] <- mean((adj - truth)^2)
        
        save.image(file=out_file)
      }
    }
  }
}
