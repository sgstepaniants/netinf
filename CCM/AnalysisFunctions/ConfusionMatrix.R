confusion_matrix <- function(mats, pred_mats) {
  nvars <- dim(mats)[1]
  num_mats <- dim(mats)[3]
  
  num_possible_mats <- 2^(nvars^2 - nvars)
  confusion_mat <- matrix(0, nrow=num_possible_mats, ncol=num_possible_mats)
  for (i in 1:num_mats) {
    mat <- mats[,,i]
    mat_index <- index_mat(mat)
    pred_mat <- pred_mats[,,i]
    pred_mat_index <- index_mat(pred_mat)
    
    confusion_mat[mat_index, pred_mat_index] <- confusion_mat[mat_index, pred_mat_index] + 1
  }
  
  return(confusion_mat)
}


index_mat <- function(mat) {
  nvars <- dim(mat)[1]
  ones_off_diag <- matrix(1, nrow=nvars, ncol=nvars)
  diag(ones_off_diag) <- 0
  power_mat <- t(matrix(cumsum(as.vector(ones_off_diag)) - 1, nvars, nvars))
  diag(power_mat) <- 0
  index <- sum(sum(mat * 2^power_mat * ones_off_diag)) + 1
  return(index)
}
