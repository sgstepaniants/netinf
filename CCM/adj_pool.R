pool <- function(mat, s, stride, type="sum") {
  # mean pool the matrix mat with a pooling filter of size s x s
  n <- dim(mat)[1];
  m <- dim(mat)[2];
  if (n < s || m < s) {
    stop('matrix dimensions are smaller than dimensions of pooling filter')
  }

  mat_pool = matrix(0, n / stride, m / stride);

  for (j in 0 : ((n - s) / stride)) {
    for (k in 0 : ((m - s) / stride)) {
      window <- mat[(stride * j + 1) : (stride * j + s), (stride * k + 1) : (stride * k + s)];
      if (type == "sum") {
        mat_pool[j + 1, k + 1] = sum(window);
      } else if (type == "mean") {
        mat_pool[j + 1, k + 1] = mean(window);
      } else {
        stop('type of matrix pool is either sum or mean')
      }
    }
  }
  
  return(mat_pool)
}
