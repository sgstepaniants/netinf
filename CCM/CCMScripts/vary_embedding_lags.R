library(rEDM)
library(lattice)
require(R.matlab)
require(nonlinearTseries)
source('ccm_helper.R')
source('CCMBaseExperiment.R')
source('AnalysisFunctions/ConfusionMatrix.R')

exp_name <- "EXPPaper1"
exp_path <- sprintf("../../KuramotoExperiments/%s", exp_name)

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
E <- 10
lib <- c(1, 100)
pred <- c(101, 250)

mat_idx <- 1
data <- data_log[,, 1, mat_idx]
mat <- true_mats[,, mat_idx]
matplot(t(data), type='l')

simplex_output <- simplex(t(data), lib, pred)
plot(simplex_output$E, simplex_output$rho, type = "l", xlab = "Embedding Dimension (E)", 
     ylab = "Forecast Skill (rho)")

num_samples <- 100
lib_sizes <- 10 : 240
xmap1 <- ccm(t(data), E=E, lib_column=1,
            target_column=2, lib_sizes=lib_sizes, tau=15,
            num_samples=num_samples, random_libs=TRUE, replace=TRUE);
xmap2 <- ccm(t(data), E=E, lib_column=2,
             target_column=1, lib_sizes=lib_sizes, tau=15,
             num_samples=num_samples, random_libs=TRUE, replace=TRUE);
xmap1_means <- ccm_means(xmap1)
xmap2_means <- ccm_means(xmap2)

plot(xmap1_means$lib_size, pmax(0, xmap1_means$rho), type = "l", col = "red", xlab = "Library Size", 
     ylab = "Cross Map Skill (rho)", ylim = c(0, 1), main = "2 causes 1")
legend(x = "topleft", legend = "1 xmap 2",
       col = c("red"), lwd = 1, bty = "n", inset = 0.02, cex = 0.8)
plot(xmap2_means$lib_size, pmax(0, xmap2_means$rho), type = "l", col = "red", xlab = "Library Size", 
     ylab = "Cross Map Skill (rho)", ylim = c(0, 1), main = "1 causes 2")
legend(x = "topleft", legend = "2 xmap 1",
       col = c("red"), lwd = 1, bty = "n", inset = 0.02, cex = 0.8)
