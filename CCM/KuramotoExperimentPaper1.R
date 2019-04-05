library(rEDM)
library(lattice)
require(R.matlab)
source('ccm_helper.R')
source('CCMBaseExperiment.R')
source('AnalysisFunctions/ConfusionMatrix.R')

exp_name <- "EXPPaper1"
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
data_log <- readMat(sprintf("%s/dataLog.mat", exp_path))[[1]]
num_mats <- dim(data_log)[4]

# Get connectivity matrices.
true_mats <- readMat(sprintf("%s/trueMats.mat", exp_path))[[1]]

# Perform CCM analysis on sample data
E <- 3
num_libs <- 10
num_trials <- 10
num_samples <- 100

exp_params <- list("E"=E, "num_libs"=num_libs, "num_trials"=num_trials, "num_samples"=num_samples)

result <- CCMBaseExperiment(data_log, true_mats, E, num_libs, num_trials, num_samples, preprocfn=cos)
pred_mats <- result$pred_mats
ccm_rho_graphs <- result$graphs

# Compute the confusion matrix
confusion_mat <- confusion_matrix(true_mats, pred_mats)
print(confusion_mat)

# Save experiment result files.
saveRDS(exp_params, sprintf("%s/exp_params.rds", result_path))
saveRDS(pred_mats, sprintf("%s/pred_mats.rds", result_path))
saveRDS(ccm_rho_graphs, sprintf("%s/ccm_rho_graphs.rds", result_path))
saveRDS(confusion_mat, sprintf("%s/confusion_mat.rds", result_path))
