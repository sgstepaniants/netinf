CCMBaseExperiment <- function(data, mats, E, num_libs, tau=1, num_trials=dim(data)[3], num_samples=100,
                              preprocfn=identity, freq=1, result_path=NULL, data_obs_idx=NULL) {
  # return: votingMats, est, tableResults
  save_data <- TRUE
  if (is.null(result_path)) {
    save_data <- FALSE
  }

  # Make directory to hold results files if one does not already exist
  if (save_data && !dir.exists(result_path)) {
    print(sprintf("Result path not found: %s", result_path))
    stop()
  }

  nvars <- dim(mats)[1] # number of variables / oscillators
  num_mats <- 1
  if (length(dim(mats)) == 3) {
    num_mats <- dim(mats)[3] # number of matrices we try
  }
  if (is.null(data_obs_idx)) {
    data_obs_idx <- matrix(TRUE, nrow=num_mats, ncol=nvars)
  }

  # Tables to hold results
  table_results_TPR <- numeric(length=num_mats)
  table_results_FPR <- numeric(length=num_mats)
  table_results_acc <- numeric(length=num_mats)

  # est holds CCM's estimate of the networks
  est <- array(0, c(nvars, nvars, num_mats))
  ccm_rho_graphs <- array(0, c(nvars, nvars, num_libs, num_trials, num_mats))

  count <- 1 # number of times we have run network inference method (so know how often to save work)

  # Loop over the networks
  for (i in 1:num_mats) {
    obs_idx <- data_obs_idx[i,]
    num_obs <- sum(obs_idx)
    
    hasMats = length(dim(mats)) == 3
    truth <- mats
    if (hasMats) {
      truth <- mats[,,i]
    }
    num_positives <- sum(truth[obs_idx, obs_idx] & !diag(nvars))
    num_negatives <- sum(!truth[obs_idx, obs_idx] & !diag(nvars))
    
    hasTrials = length(dim(data)) >= 3
    if (is.list(data)) {
      preproc_data <- preprocfn(data[[i]][[1]][obs_idx,,])
    } else if (!hasTrials && !hasMats) {
      preproc_data <- preprocfn(data[obs_idx,])
    } else if (hasTrials && !hasMats) {
      preproc_data <- preprocfn(data[obs_idx,,])
    } else {
      preproc_data <- preprocfn(data[obs_idx,,,i])
    }

    time_length <- dim(preproc_data)[2]
    delta <- floor(time_length/(num_libs + 1))
    lib_sizes <- delta * 1:num_libs

    # Run network inference on this data
    graphs <- get_ccm_rho(preproc_data, E, lib_sizes, tau, num_trials, num_samples)

    ccm_rho_graphs[,,,, i] <- graphs

    # Compute the adjacency matrix predicted by averaging the CCM trials
    voted_graphs <- apply(graphs, c(1, 2, 3), mean)
    est[,, i] <- apply(voted_graphs, c(1, 2), mean) > 0.5

    # Save results
    table_results_TPR[i] <- sum((est[,, i] + truth == 2) * !diag(nvars)) / num_positives
    table_results_FPR[i] = sum((est[,, i] - truth == 1) * !diag(nvars)) / num_negatives
    table_results_acc[i] = sum((est[,, i] == truth) * !diag(nvars)) / (num_obs^2-num_obs)

    if (save_data && count %% freq == 0) {
      # save necessary files
    }

    count <- count + 1
  }
  
  # TODO: Save whole workspace (including all those tables of results)
  if (save_data) {
    
  }
  
  table_results <- list("tpr"=table_results_TPR, "fpr"=table_results_FPR, "acc"=table_results_acc)
  result <- list("pred_mats"=est, "graphs"=ccm_rho_graphs, "table_results"=table_results)
  return(result)
}
