library(rEDM)
library(lattice)
require(R.matlab)
require(nonlinearTseries)
source('ccm_helper.R')
source('CCMBaseExperiment.R')
source('AnalysisFunctions/ConfusionMatrix.R')

exp_name <- "EXPPaper1"
exp_path <- sprintf("../../HarmonicExperiments/%s", exp_name)

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
E <- 2
lib <- c(1, 100)
pred <- c(101, 250)

data <- data_log[,, 1, 1]
matplot(t(data), type='l')

simplex_output <- simplex(t(data), lib, pred)
plot(simplex_output$E, simplex_output$rho, type = "l", xlab = "Embedding Dimension (E)", 
     ylab = "Forecast Skill (rho)")
