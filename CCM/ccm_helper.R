get_ccm_means <- function(data_func, E, names=names(data),
                          num_samples=100,
                          lib_sizes=seq(10, 80, by = 10), trials=1) {
  adj <- array(0, c(length(names), length(names), length(lib_sizes)))
  
  data <- list()
  for (t in 1:trials) {
    data[[t]] <- data_func()
    print(paste("trial", t))
  }
  
  for (i in 1:length(names)) {
    for (j in 1:length(names)) {
      sum <- numeric(length(lib_sizes))
      for (t in 1:trials) {
        xmap <- ccm(data[[t]], E = E, lib_column = names[i], 
                    target_column = names[j], lib_sizes = lib_sizes,
                    num_samples = num_samples, random_libs = TRUE, replace = TRUE);
        xmap_means <- ccm_means(xmap)
        sum <- sum + pmax(0, xmap_means$rho)
        
        #plot(xmap_means$lib_size, pmax(0, xmap_means$rho), type = "l", col = "red", xlab = "Library Size", 
        #     ylab = "Cross Map Skill (rho)", ylim = c(0, 1), main = paste(i, j))
        #legend(x = "topleft", legend = paste(i, "xmap", j),
        #       col = c("red"), lwd = 1, bty = "n", inset = 0.02, cex = 0.8)
      }
      adj[i, j,] <- sum / trials
      print(paste("xmap", i, j))
    }
  }
  
  return(adj)
}
