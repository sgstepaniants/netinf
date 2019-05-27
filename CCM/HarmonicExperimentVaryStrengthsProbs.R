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

# Load experiment simulations and ground-truth connectivity matrices
data_log <- readMat(sprintf("%s/dataLog.mat", exp_path))[[1]]
true_mats <- readMat(sprintf("%s/trueMats.mat", exp_path))[[1]]

# Perform CCM analysis on data
E <- 2
num_libs <- 1
num_samples <- 100

#window_size <- 10
#preprocfn <- function(x) {moving_average(x, window_size)}
preprocfn <- identity

# Save experiment parameters
exp_params <- list("E"=E, "num_libs"=num_libs, "num_samples"=num_samples, "preprocfn"=preprocfn)
saveRDS(exp_params, sprintf("%s/exp_params.rds", result_path))

nvars <- params$nvars
num_trials <- params$numTrials
num_mats <- params$numMats
probs <- params$probs
strengths <- params$strengths
num_probs <- params$numProbs
num_strengths <- params$numStrengths

# Register number of cores
registerDoParallel(cores=2)

# Iterate over all possible connection probabilities and spring constants
results <-
  foreach (i = 1:num_probs, .combine='cbind') %:%
    foreach (j = 1:num_strengths, .combine='cbind') %:%
      foreach (m = 1:num_mats, .combine='cbind') %dopar% {
        prob = probs[i]
        strength = strengths[j]
        
        data <- adrop(data_log[,,, i, j, m, drop=FALSE], drop=4:5)
        mat <- adrop(true_mats[,, i, j, m, drop=FALSE], drop=3:4)
        result <- CCMBaseExperiment(data, mat, E, num_libs, num_trials, num_samples, preprocfn)
      }

# Create data structures to hold experiment results
pred_mats <- array(NaN, c(nvars, nvars, num_probs, num_strengths, num_mats))
graph_log <- array(NaN, c(nvars, nvars, num_libs, num_trials, num_probs, num_strengths, num_mats))
tpr_log <- array(NaN, c(num_probs, num_strengths, num_mats))
fpr_log <- array(NaN, c(num_probs, num_strengths, num_mats))
acc_log <- array(NaN, c(num_probs, num_strengths, num_mats))
for (ind in 1:(num_probs*num_strengths*num_mats)) {
  result <- results[, ind]
  i = floor((ind-1) / (num_strengths * num_mats)) + 1
  j = (floor((ind-1) / num_mats) %% num_strengths) + 1
  m = ((ind-1) %% num_mats) + 1
  
  pred_mats[,, i, j, m] <- result$pred_mats
  graph_log[,,,, i, j, m] <- result$graphs
  table_results <- result$table_results
  tpr_log[i, j, m] <- table_results$tpr
  fpr_log[i, j, m] <- table_results$fpr
  acc_log[i, j, m] <- table_results$acc
}


# Save experiment result files.
saveRDS(pred_mats, sprintf("%s/pred_mats.rds", result_path))
writeMat(sprintf("%s/predMats.m", result_path), A = pred_mats)

saveRDS(tpr_log, sprintf("%s/tpr_log.rds", result_path))
writeMat(sprintf("%s/tprLog.m", result_path), A = tpr_log)

saveRDS(fpr_log, sprintf("%s/fpr_log.rds", result_path))
writeMat(sprintf("%s/fprLog.m", result_path), A = fpr_log)

saveRDS(acc_log, sprintf("%s/acc_log.rds", result_path))
writeMat(sprintf("%s/accLog.m", result_path), A = acc_log)

saveRDS(graph_log, sprintf("%s/graph_log.rds", result_path))
writeMat(sprintf("%s/graphLog.m", result_path), A = graph_log)


# Plot accuracy, TPR, and FPR for all connections probability and spring constant combinations
rgb.palette <- colorRampPalette(c("blue", "red"), space = "rgb")
myPanel <- function(x, y, z, ...) {
  panel.levelplot(x,y,z,...)
  panel.text(x, y, round(z,3))
}

ave_acc <- apply(acc_log, c(1, 2), function(x) {mean(x, na.rm=TRUE)})
writeMat(sprintf("%s/aveAcc.m", result_path), A = ave_acc)
levelplot(t(ave_acc), main="Average Accuracy over Simulations",
          xlab="Connection Strength", ylab="Connection Probability",
          col.regions=rgb.palette(120),
          at=seq(0, 1, length.out=120),
          panel=myPanel)

ave_tpr <- apply(tpr_log, c(1, 2), function(x) {mean(x, na.rm=TRUE)})
writeMat(sprintf("%s/aveTPR.m", result_path), A = ave_tpr)
levelplot(t(ave_tpr), main="Average TPR over Simulations",
          xlab="Connection Strength", ylab="Connection Probability",
          col.regions=rgb.palette(120),
          at=seq(0, 1, length.out=120),
          panel=myPanel)

ave_fpr <- apply(fpr_log, c(1, 2), function(x) {mean(x, na.rm=TRUE)})
writeMat(sprintf("%s/aveFPR.m", result_path), A = ave_fpr)
levelplot(t(ave_fpr), main="Average FPR over Simulations",
          xlab="Connection Strength", ylab="Connection Probability",
          col.regions=rgb.palette(120),
          at=seq(0, 1, length.out=120),
          panel=myPanel)
