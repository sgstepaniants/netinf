library(rEDM)
library(lattice)
require(R.matlab)
source('ccm_helper.R')
source('AnalysisFunctions/moving_average.R')
source('CCMBaseExperiment.R')

nvars <- 2
prob <- 0.5
spring <- 1
damping <- 0.1
pert_force <- 30

exp_name <- sprintf("EXPVaryPertsObs(nvars%d_prob%.2f_spring%.2f_damping%.2f_pertf%.2f)", nvars, prob, spring, damping, pert_force);
exp_path <- sprintf("../HarmonicExperiments/%s", exp_name);

print(exp_path)

if (!dir.exists(exp_path)) {
  print(sprintf("Data not found: %s", exp_path))
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
E <- 3
num_libs <- 10
num_trials <- 10
num_samples <- 100
num_mats <- params$numMats

window_size <- 10
preprocfn <- function(x) {moving_average(x, window_size)}

pred_mats <- array(NaN, c(nvars, nvars, num_mats, nvars, nvars))
graph_log <- array(NaN, c(nvars, nvars, num_libs, num_trials, num_mats, nvars, nvars))
norm_log <- array(NaN, c(nvars, nvars, num_mats))
tpr_log <- array(NaN, c(nvars, nvars, num_mats))
fpr_log <- array(NaN, c(nvars, nvars, num_mats))
acc_log <- array(NaN, c(nvars, nvars, num_mats))
table_results_log <- matrix(list(), nrow=nvars, ncol=nvars)

# Iterate over all possible observations and perturbations
for (num_obs in 2:nvars) {
  print(sprintf("obs %s", num_obs))
  for (num_perts in 1:nvars) {
    print(sprintf("pert %s", num_perts))

    data_log <- readMat(sprintf("%s/numobs%d/numperts%d/dataLog.mat", exp_path, num_obs, num_perts))[[1]]
    true_mats <- readMat(sprintf("%s/numobs%d/numperts%d/trueMats.mat", exp_path, num_obs, num_perts))[[1]]
    data_obs_idx <- readMat(sprintf("%s/numobs%d/numperts%d/dataObsIdx.mat", exp_path, num_obs, num_perts))[[1]] == 1
    result <- CCMBaseExperiment(data_log, true_mats, E, num_libs, num_trials, num_samples, preprocfn, data_obs_idx=data_obs_idx)
    
    pred_mats[,,, num_obs, num_perts] <- result$pred_mats
    graph_log[,,,,, num_obs, num_perts] <- result$graphs
    table_results <- result$table_results
    table_results_log[num_obs, num_perts] <- list(table_results)
    
    norm_log[num_obs, num_perts,] <- table_results$norm
    tpr_log[num_obs, num_perts,] <- table_results$tpr
    fpr_log[num_obs, num_perts,] <- table_results$fpr
    acc_log[num_obs, num_perts,] <- table_results$acc
  }
}

# Plot accuracy, TPR, and FPR for all observation and perturbation combinations
rgb.palette <- colorRampPalette(c("blue", "red"), space = "rgb")
myPanel <- function(x, y, z, ...) {
  panel.levelplot(x,y,z,...)
  panel.text(x, y, round(z,3))
}

ave_acc <- apply(acc_log, c(1, 2), function(x) {mean(x, na.rm=TRUE)})
levelplot(t(ave_acc), main="Average Accuracy over Simulations",
          xlab="Number of Perturbed Nodes", ylab="Number of Observed Nodes",
          ylim=c(nvars + 0.5, 0.5),
          col.regions=rgb.palette(120),
          at=seq(0, 1, length.out=120),
          panel=myPanel)
