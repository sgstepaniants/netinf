library(rEDM)
library(lattice)
library(abind)
library(doParallel)
require(R.matlab)
source('ccm_helper.R')
source('CCMBaseExperiment.R')
source('AnalysisFunctions/moving_average.R')

exp_name <- "VaryStrengthsProbs_Size5"
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
num_libs <- 1
num_samples <- 100
preprocfn <- identity
max_delay <- Inf
max_emb <- Inf

# Save experiment parameters
exp_params <- list("num_libs"=num_libs, "num_samples"=num_samples, "preprocfn"=preprocfn, "max_delay"=max_delay, "max_emb"=max_emb)
saveRDS(exp_params, sprintf("%s/exp_params.rds", result_path))

nvars <- params$nvars
num_trials <- params$numTrials
num_mats <- params$numMats
probs <- params$probs
num_probs <- params$numProbs
strengths <- params$strengths
num_strengths <- params$numStrengths

# Create data structures to hold experiment results
pred_mats <- array(NaN, c(nvars, nvars, num_mats, num_probs, num_strengths))
graph_log <- array(NaN, c(nvars, nvars, num_libs, num_trials, num_mats, num_probs, num_strengths))
tpr_log <- array(NaN, c(num_probs, num_strengths, num_mats))
fpr_log <- array(NaN, c(num_probs, num_strengths, num_mats))
acc_log <- array(NaN, c(num_probs, num_strengths, num_mats))

# Register number of cores
registerDoParallel(cores=12)

# Iterate over all possible connection probabilities and spring constants
results <-
  foreach (j = 1:num_probs, .combine='cbind') %:%
    foreach (k = 1:num_strengths, .combine='cbind') %:%
      foreach (m = 1:num_mats, .combine='cbind') %do% {
        print(sprintf('prob: %d, strength: %d', j, k))
        if (!file.exists(data_path)) {
          result <- list("pred_mats"=NA, "graphs"=NA, "table_results"=list("tpr"=NA, "fpr"=NA, "acc"=NA))
        } else {
          data_path <- sprintf("%s/prob%d/strength%d/mat%d/dataLog.mat", exp_path, j, k, m)
          emb_params_path <- sprintf("%s/prob%d/strength%d/mat%d/embedParams.txt", exp_path, j, k, m)
          data_log <- readMat(data_path)
          data <- data_log$noisyData
          mat <- data_log$mat
          emb_params <- embed_params(data_path, emb_params_path, max_delay, max_emb)
          result <- CCMBaseExperiment(data, mat, emb_params$E, num_libs, emb_params$tau, num_trials, num_samples, preprocfn)
        }
      }

# Create data structures to hold experiment results
pred_mats <- array(list(), c(num_probs, num_strengths, num_mats))
graph_log <- array(list(), c(num_probs, num_strengths, num_mats))
tpr_log <- array(NaN, c(num_probs, num_strengths, num_mats))
fpr_log <- array(NaN, c(num_probs, num_strengths, num_mats))
acc_log <- array(NaN, c(num_probs, num_strengths, num_mats))
for (ind in 1:(num_probs*num_strengths*num_mats)) {
  result <- results
  if (num_probs*num_strengths*num_mats > 1) {
    result <- results[, ind]
  }
  idx <- arrayInd(ind, c(num_probs, num_strengths, num_mats))
  j <- idx[1]
  k <- idx[2]
  m <- idx[3]
  
  pred_mats[j, k, m][[1]] <- result$pred_mats
  graph_log[j, k, m][[1]] <- result$graphs
  table_results <- result$table_results
  tpr_log[j, k, m] <- table_results$tpr
  fpr_log[j, k, m] <- table_results$fpr
  acc_log[j, k, m] <- table_results$acc
}


# Save experiment result files.
save(pred_mats, tpr_log, fpr_log, acc_log, graph_log, file=sprintf("%s/results.rds", result_path))
#writeMat(sprintf("%s/results.mat", result_path), A=pred_mats, B=tpr_log, C=fpr_log, D=acc_log, E=graph_log)


# Plot accuracy, TPR, and FPR for all connections probability and spring constant combinations
rgb.palette <- colorRampPalette(c("blue", "red"), space = "rgb")
myPanel <- function(x, y, z, ...) {
  panel.levelplot(x,y,z,...)
  panel.text(x, y, round(z,3))
}

ave_acc <- apply(acc_log, c(1, 2), function(x) {mean(x, na.rm=TRUE)})
levelplot(t(ave_acc), main="Average Accuracy over Simulations",
          xlab="Connection Strength", ylab="Connection Probability",
          ylim=c(num_probs + 0.5, 0.5),
          col.regions=rgb.palette(120),
          at=seq(0, 1, length.out=120),
          panel=myPanel)
