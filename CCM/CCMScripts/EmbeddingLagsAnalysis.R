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

exp_name <- "VaryPertsObs"
exp_path <- sprintf("../HarmonicExperiments/EXP%s", exp_name)

print(exp_path)

if (!dir.exists(exp_path)) {
  print(sprintf("Data not found: %s", exp_path))
  stop()
}

max_delay = 100;
max_emb = 20;

# Read data simulation parameters
filepath <- sprintf("%s/dataLog.mat", exp_path)
data <- readMat(filepath)[[1]]
code = c("addpath('../DataScripts/SimulateData')",
         "addpath('../mdembedding')",
         sprintf("maxDelay = %d; maxEmb = %d;", max_delay, max_emb),
         sprintf("load('%s');", filepath),
         "nobs = size(noisyData, 2);",
         "tau = round(mdDelay(noisyData.', 'maxLag', maxDelay, 'plottype', 'none'));",
         "currMaxEmb = min(floor(nobs / tau), maxEmb);",
         "[fnnPercent, Es] = mdFnn(noisyData(1, :).', tau, 'maxEmb', currMaxEmb, 'doPlot', 1);",
         "E = findElbow(Es, fnnPercent);",
         "save('embedParams.txt', 'tau', 'E', '-ascii')")
res = run_matlab_code(code)
