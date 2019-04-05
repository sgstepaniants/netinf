CCMBaseExperiment <- function(data, mats, E, num_libs, num_trials=dim(data)[3], num_samples=100,
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
  num_mats <- dim(mats)[3] # number of matrices we try
  if (is.null(data_obs_idx)) {
    data_obs_idx <- matrix(TRUE, nrow=num_mats, ncol=nvars)
  }

  # Tables to hold results
  table_results_norm <- numeric(length=num_mats)
  table_results_TPR <- numeric(length=num_mats)
  table_results_FPR <- numeric(length=num_mats)
  table_results_acc <- numeric(length=num_mats)
  
  #table_results_norm <- matrix(0, nrow=num_trials, ncol=num_mats)
  #table_results_norm_voting <- numeric(length=num_mats)
  #table_results_TPR <- matrix(0, nrow=num_trials, ncol=num_mats)
  #table_results_TPR_voting <- numeric(length=num_mats)
  #table_results_FPR <- matrix(0, nrow=num_trials, ncol=num_mats)
  #table_results_FPR_voting <- numeric(length=num_mats)
  #table_results_acc <- matrix(0, nrow=num_trials, ncol=num_mats)
  #table_results_acc_voting <- numeric(length=num_mats)

  #voting_mats <- array(0, c(nvars, nvars, num_mats))

  # est holds CCM's estimate of the networks
  est <- array(0, c(nvars, nvars, num_mats))
  ccm_rho_graphs <- array(0, c(nvars, nvars, num_libs, num_trials, num_mats))

  count <- 1 # number of times we have run network inference method (so know how often to save work)

  # Loop over the networks
  for (i in 1:num_mats) {
    obs_idx <- data_obs_idx[i,]
    num_obs <- sum(obs_idx)

    truth <- mats[,,i]
    num_positives <- sum((truth[obs_idx, obs_idx] * !diag(nvars)))
    num_negatives <- sum(!truth[obs_idx, obs_idx] * !diag(nvars))

    if (is.list(data)) {
      preproc_data <- preprocfn(data[[i]][[1]][obs_idx,,])
    } else {
      preproc_data <- preprocfn(data[obs_idx,,,i])
    }

    time_length <- dim(preproc_data)[2]
    delta <- floor(time_length/(num_libs+1))
    lib_sizes <- seq(delta, time_length-delta, by=delta)

    # Run network inference on this data
    ptm <- proc.time()
    graphs <- get_ccm_rho(preproc_data, E, lib_sizes, num_trials, num_samples)
    print(proc.time() - ptm)

    ccm_rho_graphs[,,,, i] <- graphs

    # Compute the adjacency matrix predicted by averaging the CCM trials
    voted_graphs <- apply(graphs, c(1, 2, 4), mean)
    est[,, i] <- apply(voted_graphs, c(1, 2), mean) > 0.5

    # Save results
    table_results_norm[i] <- norm(est[,, i] - truth) / norm(truth)
    table_results_TPR[i] <- sum((est[,, i] + truth == 2) * !diag(nvars)) / num_positives
    table_results_FPR[i] = sum((est[,, i] - truth == 1) * !diag(nvars)) / num_negatives
    table_results_acc[i] = sum((est[,, i] == truth) * !diag(nvars)) / (num_obs^2-num_obs)

    #tableResultsNormVoting(j) = norm(votingMat - truth) / norm(truth);
    #tableResultsTPRVoting(j) = sum((votingMat + truth == 2) * ~eye(nvars)) / num_positives
    #tableResultsFPRVoting(j) = sum((votingMat - truth == 1) * ~eye(nvars)) / num_negatives
    #tableResultsAccVoting(j) = sum((votingMat == truth) * ~eye(nvars)) / (numObs^2-numObs)

    if (save_data && count %% freq == 0) {
      # save necessary files
    }

    count <- count + 1
  }

  #tableResults = struct('norm', tableResultsNorm, 'normVoting', tableResultsNormVoting, ...
  #'tpr', tableResultsTPR, 'tprVoting', tableResultsTPRVoting, ...
  #'fpr', tableResultsFPR, 'fprVoting', tableResultsFPRVoting, ...
  #'acc', tableResultsAcc, 'accVoting', tableResultsAccVoting, ...
  #'diagnostics', tableResultsDiagnostics);
  
  # TODO: Save whole workspace (including all those tables of results)
  if (save_data) {
    
  }
  
  table_results <- list("norm"=table_results_norm, "tpr"=table_results_TPR, "fpr"=table_results_FPR, "acc"=table_results_acc)
  result <- list("pred_mats"=est, "graphs"=ccm_rho_graphs, "table_results"=table_results)
  return(result)
}
