library(rEDM)
library(lattice)
require(R.matlab)
library(doParallel)
source('ccm_helper.R')
source('CCMBaseExperiment.R')
source('AnalysisFunctions/ConfusionMatrix.R')

exp_name <- "EXPConfusionMatrix_Size2"
exp_path <- sprintf("../HarmonicExperiments/%s", exp_name)

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
nvars <- 2
data_log <- readMat(sprintf("%s/dataLog.mat", exp_path))[[1]]
num_mats <- dim(data_log)[4]

# Get connectivity matrices.
true_mats <- readMat(sprintf("%s/trueMats.mat", exp_path))[[1]]

# Perform CCM analysis on sample data
Es <- seq(2, 15, 2)
taus <- c(1, 5, 10, 15)
num_Es <- length(Es)
num_taus <- length(taus)
num_libs <- 10
num_trials <- 1
num_samples <- 100

exp_params <- list("Es"=Es, "num_libs"=num_libs, "taus"=taus, "num_trials"=num_trials, "num_samples"=num_samples)
saveRDS(exp_params, sprintf("%s/exp_params.rds", result_path))

# Register number of cores
registerDoParallel(cores=12)

results <-
  foreach (i = 1 : num_Es, .combine='cbind') %:%
    foreach (j = 1 : num_taus, .combine='cbind') %dopar% {
      print(sprintf("num_E: %d, num_tau: %d", i, j))
      E <- Es[i];
      tau <- taus[j]
      
      result <- CCMBaseExperiment(data_log, true_mats, E, num_libs, tau, num_trials, num_samples)
    }

graph_log <- array(NaN, c(nvars, nvars, num_libs, num_trials, num_mats, num_Es, num_taus))
pred_mats <- array(NaN, c(nvars, nvars, num_mats, num_Es, num_taus))
tpr_log <- array(NaN, c(num_Es, num_taus, num_mats))
fpr_log <- array(NaN, c(num_Es, num_taus, num_mats))
acc_log <- array(NaN, c(num_Es, num_taus, num_mats))
for (ind in 1:(num_Es*num_taus)) {
  result <- results[, ind]
  idx <- arrayInd(ind, c(num_Es, num_taus))
  i <- idx[1]
  j <- idx[2]
  
  pred_mats[,,, i, j] <- result$pred_mats
  graph_log[,,,,, i, j] <- result$graphs
  table_results <- result$table_results
  tpr_log[i, j,] <- table_results$tpr
  fpr_log[i, j,] <- table_results$fpr
  acc_log[i, j,] <- table_results$acc
}

ave_accs <- apply(acc_log, c(1, 2), mean)

# Compute the confusion matrix
confusion_mat <- confusion_matrix(true_mats, pred_mats[,,, 1, 4])
print(confusion_mat)

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


# Plot directional causality rho graphs
ave_rho_graphs <- apply(graph_log, c(1, 2, 3, 5, 6, 7), mean)

maxInd <- 100;
rho_graphs_none <- apply(ave_rho_graphs[,,,1:maxInd,,,drop=FALSE], c(1, 2, 3, 5, 6), mean)
rho_graphs_1causes2 <- apply(ave_rho_graphs[,,,101:(maxInd+100),,,drop=FALSE], c(1, 2, 3, 5, 6), mean)
rho_graphs_2causes1 <- apply(ave_rho_graphs[,,,201:(maxInd+200),,,drop=FALSE], c(1, 2, 3, 5, 6), mean)
rho_graphs_both <- apply(ave_rho_graphs[,,,301:(maxInd+300),,,drop=FALSE], c(1, 2, 3, 5, 6), mean)

ind_E <- 1
ind_tau <- 4
delta <- floor(250/(num_libs + 1))
lib_sizes <- delta * 1:num_libs

# No causality
plot(lib_sizes, rho_graphs_none[1, 2,, ind_E, ind_tau], type = "l", col = "red", xlab = "Library Size", 
     ylab = "Cross Map Skill (rho)", ylim = c(0, 1), main = "No causality")
lines(lib_sizes, rho_graphs_none[2, 1,, ind_E, ind_tau], type = "l", col = "blue", xlab = "Library Size", 
      ylab = "Cross Map Skill (rho)", ylim = c(0, 1))
legend(x = "topleft", col = c("red", "blue"), legend = c("1 xmap 2", "2 xmap 1"),
       lwd = 2, inset = 0.02, bty = "n", cex = 0.8)

# 1 causes 2
plot(lib_sizes, rho_graphs_1causes2[1, 2,, ind_E, ind_tau], type = "l", col = "red", xlab = "Library Size", 
     ylab = "Cross Map Skill (rho)", ylim = c(0, 1), main = "1 -> 2")
lines(lib_sizes, rho_graphs_1causes2[2, 1,, ind_E, ind_tau], type = "l", col = "blue", xlab = "Library Size", 
      ylab = "Cross Map Skill (rho)", ylim = c(0, 1))
legend(x = "topleft", col = c("red", "blue"), legend = c("1 xmap 2", "2 xmap 1"),
       lwd = 2, inset = 0.02, bty = "n", cex = 0.8)

# 2 causes 1
plot(lib_sizes, rho_graphs_2causes1[1, 2,, ind_E, ind_tau], type = "l", col = "red", xlab = "Library Size", 
     ylab = "Cross Map Skill (rho)", ylim = c(0, 1), main = "1 <- 2")
lines(lib_sizes, rho_graphs_2causes1[2, 1,, ind_E, ind_tau], type = "l", col = "blue", xlab = "Library Size", 
      ylab = "Cross Map Skill (rho)", ylim = c(0, 1))
legend(x = "topleft", col = c("red", "blue"), legend = c("1 xmap 2", "2 xmap 1"),
       lwd = 2, inset = 0.02, bty = "n", cex = 0.8)

# Bidirectional causality
plot(lib_sizes, rho_graphs_both[1, 2,, ind_E, ind_tau], type = "l", col = "red", xlab = "Library Size", 
     ylab = "Cross Map Skill (rho)", ylim = c(0, 1), main = "1 <-> 2")
lines(lib_sizes, rho_graphs_both[2, 1,, ind_E, ind_tau], type = "l", col = "blue", xlab = "Library Size", 
     ylab = "Cross Map Skill (rho)", ylim = c(0, 1))
legend(x = "topleft", col = c("red", "blue"), legend = c("1 xmap 2", "2 xmap 1"),
       lwd = 2, inset = 0.02, bty = "n", cex = 0.8)
