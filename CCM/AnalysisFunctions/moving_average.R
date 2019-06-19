require(zoo)

moving_average <- function(data, window_size) {
  func <- function(x) {return(rollapply(x, window_size, mean, by=window_size))}
  if (length(dim(data)) < 3) {
    filtered_data <- t(apply(data, 1, func))
  } else {
    filtered_data <- aperm(apply(data, c(1, 3), func), c(2, 1, 3))
  }
  return(filtered_data)
}
