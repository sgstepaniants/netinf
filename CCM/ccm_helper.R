get_ccm_rho <- function(data_func, E,
                          names=names(data),
                          num_samples=100,
                          lib_sizes=seq(10, 80, by = 10),
                          trials=1) {
  ccm_rho_graphs <- array(0, c(length(names), length(names), length(lib_sizes), trials))
  
  data <- list()
  for (t in 1:trials) {
    data[[t]] <- data_func()
    print(paste("trial", t))
  }
  
  for (i in 1:length(names)) {
    for (j in 1:length(names)) {
      sum_rho_graph <- numeric(length(lib_sizes))
      for (t in 1:trials) {
        xmap <- ccm(data[[t]], E = E, lib_column = names[i], 
                    target_column = names[j], lib_sizes = lib_sizes,
                    num_samples = num_samples, random_libs = TRUE, replace = TRUE);
        xmap_means <- ccm_means(xmap)
        ccm_rho_graphs[i, j,, t] <- pmax(0, xmap_means$rho)
        
        #plot(xmap_means$lib_size, pmax(0, xmap_means$rho), type = "l", col = "red", xlab = "Library Size", 
        #     ylab = "Cross Map Skill (rho)", ylim = c(0, 1), main = paste(i, j))
        #legend(x = "topleft", legend = paste(i, "xmap", j),
        #       col = c("red"), lwd = 1, bty = "n", inset = 0.02, cex = 0.8)
      }
      print(paste("xmap", i, j))
    }
  }
  
  return(ccm_rho_graphs)
}


get_ccm_rho_hankel <- function(data_func, E, rank,
                          names=names(data),
                          num_samples=100,
                          lib_sizes=seq(10, 80, by = 10),
                          trials=1) {
  ccm_rho_graphs <- array(0, c(length(names), length(names), length(lib_sizes), trials))
  
  # store the svd of all the Hankel matrices per trial
  svd_data <- list()
  for (t in 1:trials) {
    data_hankel <- data_func()
    # the Hankel matrix is generally too large to process so take the SVD
    s <- svd(data_hankel, nu = rank, nv = rank);
    svd_data[[t]] <- s;
    print(paste("trial", t))
  }
  
  for (i in 1:length(names)) {
    for (j in 1:length(names)) {
      sum <- numeric(length(lib_sizes))
      for (t in 1:trials) {
        # get the U eigenvalue time series of the Hankel matrix
        # for the current trial
        u_data <- cbind(seq(1, tmax, by=deltat), s$u);
        colnames(u_data) <- c("time", names);
        xmap <- ccm(u_data, E = E, lib_column = names[i], 
                    target_column = names[j], lib_sizes = lib_sizes,
                    num_samples = num_samples, random_libs = TRUE, replace = TRUE);
        xmap_means <- ccm_means(xmap)
        ccm_rho_graphs[i, j,, t] <- pmax(0, xmap_means$rho)
        
        plot(xmap_means$lib_size, pmax(0, xmap_means$rho), type = "l", col = "red", xlab = "Library Size", 
             ylab = "Cross Map Skill (rho)", ylim = c(0, 1), main = paste(i, j))
        legend(x = "topleft", legend = paste(i, "xmap", j),
               col = c("red"), lwd = 1, bty = "n", inset = 0.02, cex = 0.8)
      }
      print(paste("xmap", i, j))
    }
  }
  
  return(list("ccm_rho_graphs" = ccm_rho_graphs, "svd_data" = svd_data))
}
