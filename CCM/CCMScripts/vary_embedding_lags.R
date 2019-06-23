library(rEDM)
library(lattice)
require(R.matlab)
require(nonlinearTseries)
setwd("~/netinf/CCM/CCMScripts")
source('../ccm_helper.R')
source('../CCMBaseExperiment.R')
source('../AnalysisFunctions/ConfusionMatrix.R')
source('../AnalysisFunctions/moving_average.R')

exp_name <- "EXPPertVarySizeForcingStrengths"
exp_path <- sprintf("../../HarmonicExperiments/%s", exp_name)

print(exp_path)

if (!dir.exists(exp_path)) {
  print(sprintf("Data not found: %s", exp_path))
  stop()
}

window_size <- 10
preprocfn <- function(x) {moving_average(x, window_size)}

# Get trajectory data.
j <- 1
k <- 1
l <- 1
m <- 1
data_log <- readMat(sprintf("%s/size%d/force%d/strength%d/mat%d/dataLog.mat", exp_path, j, k, l, m))
data <- preprocfn(data_log$noisyData)
nobs <- dim(data)[2]
mat <- data_log$mat

matplot(t(data), type='l')

# Perform CCM analysis on sample data
E <- 10
tau <- 1

lib <- c(1, round(nobs / 3))
pred <- c(round(nobs / 3) + 1, nobs)

simplex_output <- simplex(t(data), lib, pred)
plot(simplex_output$E, simplex_output$rho, type = "l", xlab = "Embedding Dimension (E)", 
     ylab = "Forecast Skill (rho)")

num_samples <- 100
lib_sizes <- seq(10, nobs, 10)

ind1 <- 1;
ind2 <- 2;
xmap1 <- ccm(t(data), E=E, lib_column=ind1,
            target_column=ind2, lib_sizes=lib_sizes, tau=tau,
            num_samples=num_samples, random_libs=TRUE, replace=TRUE);
xmap2 <- ccm(t(data), E=E, lib_column=ind2,
             target_column=ind1, lib_sizes=lib_sizes, tau=tau,
             num_samples=num_samples, random_libs=TRUE, replace=TRUE);
xmap1_means <- ccm_means(xmap1)
xmap2_means <- ccm_means(xmap2)

plot(xmap1_means$lib_size, pmax(0, xmap1_means$rho), type = "l", col = "red", xlab = "Library Size", 
     ylab = "Cross Map Skill (rho)", ylim = c(0, 1), main = paste(ind2, "causes", ind1))
legend(x = "topleft", legend = paste(ind1, "xmap", ind2),
       col = c("red"), lwd = 1, bty = "n", inset = 0.02, cex = 0.8)
plot(xmap2_means$lib_size, pmax(0, xmap2_means$rho), type = "l", col = "red", xlab = "Library Size", 
     ylab = "Cross Map Skill (rho)", ylim = c(0, 1), main = paste(ind1, "causes", ind2))
legend(x = "topleft", legend = paste(ind2, "xmap", ind1),
       col = c("red"), lwd = 1, bty = "n", inset = 0.02, cex = 0.8)
