library(rEDM)
library(lattice)
library(abind)
library(doParallel)
require(R.matlab)
setwd("~/netinf/CCM")
source('ccm_helper.R')
source('CCMBaseExperiment.R')
source('AnalysisFunctions/moving_average.R')

exp_name <- "PertVarySizeForcingStrengths"
exp_path <- sprintf("../HarmonicExperiments/EXP%s", exp_name)

print(exp_path)

if (!dir.exists(exp_path)) {
  print(sprintf("Data not found: %s", exp_path))
  stop()
}

result_path <- sprintf("%s/CCMResults", exp_path)
if (!dir.exists(result_path)) {
  dir.create(result_path)
} else {
  m <- readline(prompt=sprintf("%s\n already exists, would you like to continue and overwrite these results (Y/N): ", result_path))
  if (toupper(m) == "N") {
    stop()
  }
}


# Read data simulation parameters
params <- readMat(sprintf("%s/params.mat", exp_path))

# Perform CCM analysis on data
E <- 2
num_libs <- 1
num_trials <- 1
num_samples <- 100

window_size <- 10
preprocfn <- function(x) {moving_average(x, window_size)}
#preprocfn <- identity

# Save experiment parameters
exp_params <- list("E"=E, "num_libs"=num_libs, "num_samples"=num_samples, "preprocfn"=preprocfn)
saveRDS(exp_params, sprintf("%s/exp_params.rds", result_path))

num_sizes <- params$numSizes
num_forces <- params$numForces
num_strengths <- params$numStrengths
num_mats <- params$numMats

# Register number of cores
registerDoParallel(cores=2)

max_E <- 10

# Iterate over all possible connection probabilities and spring constants
count <- 1
results <-
  foreach (j = num_sizes:1, .combine='cbind') %:%
    foreach (k = 1:num_forces, .combine='cbind') %:%
      foreach (l = 1:num_strengths, .combine='cbind') %:%
        foreach (m = 1:num_mats, .combine='cbind') %:%
          foreach (E = 1:max_E, .combine='cbind') %do% {
            print(count)
            data <- readMat(sprintf("%s/size%d/force%d/dataLog.mat", exp_path, j, k))$currDataLog[[(m - 1) * num_strengths + l]][[1]]
            time_length <- dim(data)[2]
            lib <- c(1, floor(time_length / 3))
            pred <- c(floor(time_length / 3) + 1, time_length)
            result <- simplex(t(data), lib, pred)
            count <- count + 1
          }

graph_log <- array(list(), c(num_sizes, num_forces, num_strengths, num_mats))
optimal_E <- array(NaN, c(num_sizes, num_forces, num_strengths, num_mats))
for (ind in 1:(num_sizes*num_forces*num_strengths*num_mats)) {
  result <- results
  if (num_sizes*num_forces*num_strengths*num_mats > 1) {
    result <- results[, ind]
  }
  idx <- arrayInd(ind, c(num_sizes, num_forces, num_strengths, num_mats))
  j <- idx[1]
  k <- idx[2]
  l <- idx[3]
  m <- idx[4]
  
  graph_log[j, k, l, m][[1]] <- result$rho
  optimal_E[j, k, l, m] <- which(result == max(result$rho))
}
