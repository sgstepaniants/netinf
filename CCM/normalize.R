normalize <- function(mat) {
  return((mat - min(mat)) / (max(mat)-min(mat)))
}
