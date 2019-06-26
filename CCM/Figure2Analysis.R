library(rEDM)
library(lattice)
require(R.matlab)
source('ccm_helper.R')
source('CCMBaseExperiment.R')
source('AnalysisFunctions/ConfusionMatrix.R')

exp_name <- "EXPNetworkInferenceAccuracy"
exp_path <- sprintf("../KuramotoExperiments/%s", exp_name)

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


# Get trajectory data.
data_log <- readMat(sprintf("%s/dataLog.mat", exp_path))
data1 <- data_log[[1]]
data2 <- data_log[[2]]
data3 <- data_log[[3]]

# Get connectivity matrix.
true_mat <- readMat(sprintf("%s/trueMat.mat", exp_path))[[1]]

# Perform CCM analysis on sample data
E <- 5
tau <- 5
num_libs <- 10
num_trials <- 100
num_samples <- 100

exp_params <- list("E"=E, "num_libs"=num_libs, "tau"=tau, "num_trials"=num_trials, "num_samples"=num_samples)
saveRDS(exp_params, sprintf("%s/exp_params.rds", result_path))

result1 <- CCMBaseExperiment(data1, true_mat, E, num_libs, tau, num_trials, num_samples)
pred_mat1 <- drop(result1$pred_mats)
graph1 <- drop(result1$graph)
table_results1 <- result1$table_results
tpr1 <- table_results1$tpr
fpr1 <- table_results1$fpr
acc1 <- table_results1$acc

# Save experiment result files.
saveRDS(pred_mats, sprintf("%s/pred_mats.rds", result_path))
saveRDS(ccm_rho_graphs, sprintf("%s/ccm_rho_graphs.rds", result_path))
saveRDS(confusion_mat, sprintf("%s/confusion_mat.rds", result_path))
