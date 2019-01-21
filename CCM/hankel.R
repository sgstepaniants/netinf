hankel <- function(data, stacks) {
  # remove the first column of the data matrix which
  # represents time
  data <- data[, -1];
  
  # stack data matrix into Hankel time series matrix
  n <- ncol(data);
  m <- nrow(data);
  
  data_hankel = matrix(0, nrow = m - stacks + 1, ncol = n * stacks);
  
  for (k in 1 : (m - stacks + 1)) {
    data_slice = data[k : (k + stacks - 1),];
    data_hankel[k,] = c(data_slice);
  }
  return(data_hankel);
}
