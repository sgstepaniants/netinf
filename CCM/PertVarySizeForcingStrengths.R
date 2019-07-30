library(rEDM)
library(lattice)
library(abind)
library(doParallel)
require(R.matlab)
require(matlabr)
setwd("~/netinf/CCM")
source('ccm_helper.R')
source('CCMBaseExperiment.R')
source('AnalysisFunctions/moving_average.R')

exp_name <- "PertVarySizeForcingStrengths"
exp_path <- sprintf("../KuramotoExperiments/EXP%s", exp_name)

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

num_trials <- 1 #params$numTrials
num_mats <- params$numMats
network_sizes <- params$networkSizes
num_sizes <- params$numSizes
forces <- params$forces
num_forces <- params$numForces
strengths <- params$strengths
num_strengths <- params$numStrengths

# Register number of cores
registerDoParallel(cores=25)

# Iterate over all possible connection probabilities and spring constants
Es <- array(NaN, c(num_sizes, num_forces, num_strengths, num_mats))
taus <- array(NaN, c(num_sizes, num_forces, num_strengths, num_mats))
results <-
  foreach (j = 1:num_sizes, .combine='cbind') %:%
    foreach (k = 1:num_forces, .combine='cbind') %:%
      foreach (l = 1:num_strengths, .combine='cbind') %:%
        foreach (m = 1:num_mats, .combine='cbind') %dopar% {
          print(sprintf('size: %d, force: %d, strength: %d, mat: %d', j, k, l, m))
          data_path <- sprintf("%s/size%d/force%d/strength%d/mat%d/dataLog.mat", exp_path, j, k, l, m)
          emb_params_path <- sprintf("%s/size%d/force%d/strength%d/mat%d/embedParams.txt", exp_path, j, k, l, m)
          data_log <- readMat(data_path)
          data <- data_log$noisyData
          mat <- data_log$mat
          emb_params <- embed_params(data_path, emb_params_path, max_delay, max_emb)
          Es[j, k, l, m] <- emb_params$E
          taus[j, k, l, m] <- emb_params$tau
          
          result <- CCMBaseExperiment(data, mat, emb_params$E, num_libs, emb_params$tau, num_trials, num_samples, preprocfn)
          save(result, file=sprintf("%s/size%d/force%d/strength%d/mat%d/result.rds", exp_path, j, k, l, m))
          result
        }

# Create data structures to hold experiment results
pred_mats <- array(list(), c(num_sizes, num_forces, num_strengths, num_mats))
graph_log <- array(list(), c(num_sizes, num_forces, num_strengths, num_mats))
tpr_log <- array(NaN, c(num_sizes, num_forces, num_strengths, num_mats))
fpr_log <- array(NaN, c(num_sizes, num_forces, num_strengths, num_mats))
acc_log <- array(NaN, c(num_sizes, num_forces, num_strengths, num_mats))
for (ind in 1:(num_sizes*num_forces*num_strengths*num_mats)) {
  result <- results
  if (num_sizes*num_forces*num_strengths*num_mats > 1) {
    result <- results[, ind]
  }
  idx <- arrayInd(ind, c(num_mats, num_strengths, num_forces, num_sizes))
  m <- idx[1]
  l <- idx[2]
  k <- idx[3]
  j <- idx[4]
  
  pred_mats[j, k, l, m][[1]] <- result$pred_mats
  graph_log[j, k, l, m][[1]] <- result$graphs
  table_results <- result$table_results
  tpr_log[j, k, l, m] <- table_results$tpr
  fpr_log[j, k, l, m] <- table_results$fpr
  acc_log[j, k, l, m] <- table_results$acc
}


# Save experiment result files.
save(pred_mats, tpr_log, fpr_log, acc_log, Es, taus, graph_log, file=sprintf("%s/results.rds", result_path))
writeMat(sprintf("%s/results.mat", result_path), tprLog=tpr_log, fprLog=fpr_log, accLog=acc_log)


# Plot accuracy, TPR, and FPR for all connections probability and spring constant combinations
rgb.palette <- colorRampPalette(c("blue", "red"), space = "rgb")
myPanel <- function(x, y, z, ...) {
  panel.levelplot(x,y,z,...)
  panel.text(x, y, round(z,3))
}

forceInd <- 1

ave_acc <- apply(acc_log[, forceInd,,], c(1, 2), function(x) {mean(x, na.rm=TRUE)})
writeMat(sprintf("%s/aveAcc.m", result_path), A = ave_acc)
levelplot(t(ave_acc), main="Average Accuracy over Simulations",
          xlab="Connection Strength", ylab="Connection Probability",
          col.regions=rgb.palette(120),
          at=seq(0, 1, length.out=120),
          panel=myPanel)

ave_tpr <- apply(tpr_log[, forceInd,,], c(1, 2), function(x) {mean(x, na.rm=TRUE)})
writeMat(sprintf("%s/aveTPR.m", result_path), A = ave_tpr)
levelplot(t(ave_tpr), main="Average TPR over Simulations",
          xlab="Connection Strength", ylab="Connection Probability",
          col.regions=rgb.palette(120),
          at=seq(0, 1, length.out=120),
          panel=myPanel)

ave_fpr <- apply(fpr_log[, forceInd,,], c(1, 2), function(x) {mean(x, na.rm=TRUE)})
writeMat(sprintf("%s/aveFPR.m", result_path), A = ave_fpr)
levelplot(t(ave_fpr), main="Average FPR over Simulations",
          xlab="Connection Strength", ylab="Connection Probability",
          col.regions=rgb.palette(120),
          at=seq(0, 1, length.out=120),
          panel=myPanel)
