library(rEDM)
library(lattice)
library(abind)
require(R.matlab)
source('ccm_helper.R')
source('CCMBaseExperiment.R')
source('AnalysisFunctions/moving_average.R')

exp_name <- "Paper2"
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

# Load experiment simulations and ground-truth connectivity matrices
data_log <- readMat(sprintf("%s/dataLog.mat", exp_path))[[1]]
true_mats <- readMat(sprintf("%s/trueMats.mat", exp_path))[[1]]

# Perform CCM analysis on data
E <- 2
num_libs <- 1
num_samples <- 100

window_size <- 10
#preprocfn <- function(x) {moving_average(x, window_size)}
preprocfn <- identity

# Save experiment parameters
exp_params <- list("E"=E, "num_libs"=num_libs, "num_trials"=num_trials,
                   "num_samples"=num_samples, "preprocfn"=preprocfn, "window_size"=window_size)
saveRDS(exp_params, sprintf("%s/exp_params.rds", result_path))

nvars <- params$nvars
num_trials <- 100 #params$numTrials
num_mats <- params$numMats
probs <- params$probs
strengths <- params$strengths
num_probs <- params$numProbs
num_strengths <- params$numStrengths

# Create data structures to hold experiment results
pred_mats <- array(NaN, c(nvars, nvars, num_mats, num_probs, num_strengths))
graph_log <- array(NaN, c(nvars, nvars, num_libs, num_trials, num_mats, num_probs, num_strengths))
tpr_log <- array(NaN, c(num_probs, num_strengths, num_mats))
fpr_log <- array(NaN, c(num_probs, num_strengths, num_mats))
acc_log <- array(NaN, c(num_probs, num_strengths, num_mats))
table_results_log <- matrix(list(), nrow=num_probs, ncol=num_strengths)

# Iterate over all possible connection probabilities and spring constants
for (i in 1:num_probs) {
  prob = probs[i]
  print(sprintf("prob %f", prob))
  for (j in 1:num_strengths) {
    strength = strengths[j]
    print(sprintf("strength %f", strength))
    
    data <- adrop(data_log[,,,, i, j, drop=FALSE], drop=5:6)
    mats <- adrop(true_mats[,,, i, j, drop=FALSE], drop=4:5)
    result <- CCMBaseExperiment(data, mats, E, num_libs, num_trials, num_samples, preprocfn)
    
    pred_mats[,,, i, j] <- result$pred_mats
    graph_log[,,,,, i, j] <- result$graphs
    table_results <- result$table_results
    table_results_log[i, j] <- list(table_results)
    
    tpr_log[i, j,] <- table_results$tpr
    fpr_log[i, j,] <- table_results$fpr
    acc_log[i, j,] <- table_results$acc
  }
}

# Save experiment result files.
saveRDS(pred_mats, sprintf("%s/pred_mats.rds", result_path))
saveRDS(graph_log, sprintf("%s/graph_log.rds", result_path))
saveRDS(table_results_log, sprintf("%s/table_results_log.rds", result_path))
saveRDS(tpr_log, sprintf("%s/tpr_log.rds", result_path))
saveRDS(fpr_log, sprintf("%s/fpr_log.rds", result_path))
saveRDS(acc_log, sprintf("%s/acc_log.rds", result_path))

# Plot accuracy, TPR, and FPR for all connections probability and spring constant combinations
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
