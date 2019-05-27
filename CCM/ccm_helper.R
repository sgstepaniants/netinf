# TODO: Maybe add a parameter for coordinate randomization.
get_ccm_rho <- function(data, E, lib_sizes, num_trials=dim(data)[3], num_samples=100) {
  nvars <- dim(data)[1]
  ccm_rho_graphs <- array(0, c(nvars, nvars, length(lib_sizes), num_trials))
  
  for (i in 1:nvars) {
    for (j in 1:nvars) {
      if (i != j) {
        sum_rho_graph <- numeric(length(lib_sizes))
        for (t in 1:num_trials) {
          xmap <- ccm(t(data[,,t]), E=E, lib_column=i, 
                      target_column=j, lib_sizes=lib_sizes,
                      num_samples=num_samples, random_libs=TRUE, replace=TRUE);
          xmap_means <- ccm_means(xmap)
          
          ccm_rho_graphs[i, j,, t] <- pmax(0, xmap_means$rho)
          
          #plot(xmap_means$lib_size, pmax(0, xmap_means$rho), type = "l", col = "red", xlab = "Library Size", 
          #     ylab = "Cross Map Skill (rho)", ylim = c(0, 1), main = paste(j, "causes", i))
          #legend(x = "topleft", legend = paste(i, "xmap", j),
          #       col = c("red"), lwd = 1, bty = "n", inset = 0.02, cex = 0.8)
        }
        #print(paste(j, "causes", i))
      }
    }
  }
  
  return(ccm_rho_graphs)
}


# TODO: Make get_adj returns two adjacency matrices, one containing complexity estimates
# and the other containing directionality estimates.
#get_adj <- function(ccm_rho_graphs) {
#  # ccm_rho_graphs is a 4 dimensional array where the first and second
#  # dimensions represent the nodes in the system, the third
#  # dimension is the number of library sizes tried, and the fourth
#  # dimension is the number of simulation trials ran.
#  n <- dim(ccm_rho_graphs)[1]
#  num_trials <- dim(ccm_rho_graphs)[4]
#  
#  adj_opt_rho <- array(0, dim = c(n, n, num_trials))
#  max_rho <- apply(ccm_rho_graphs, c(1, 2, 4), max)
#  
#  # Search through rho values and find the one that represents 95%
#  # of the correlation maximum for each graph.
#  for (i in 1 : n) {
#    for (j in 1 : n) {
#      for (t in 1 : num_trials) {
#        idx <- which(ccm_rho_graphs[i, j,, t] >= 0.95 * max_rho[i, j, t])[1]
#        adj_opt_rho[i, j, t] <- ccm_rho_graphs[i, j, idx, t]
#      }
#    }
#  }
#  
#  return(adj_opt_rho)
#}


# get_ccm_rho_hankel <- function(data_func, E, rank,
#                                names=names(data),
#                                num_samples=100,
#                                lib_sizes=seq(10, 80, by = 10),
#                                num_trials=1) {
#   ccm_rho_graphs <- array(0, c(length(names), length(names), length(lib_sizes), num_trials))
#   
#   # store the svd of all the Hankel matrices per trial
#   svd_data <- list()
#   for (t in 1:num_trials) {
#     data_hankel <- data_func()
#     # the Hankel matrix is generally too large to process so take the SVD
#     s <- svd(data_hankel, nu = rank, nv = rank);
#     svd_data[[t]] <- s;
#     print(paste("trial", t))
#   }
#   
#   for (i in 1:length(names)) {
#     for (j in 1:length(names)) {
#       sum <- numeric(length(lib_sizes))
#       for (t in 1:num_trials) {
#         # get the U eigenvalue time series of the Hankel matrix
#         # for the current trial
#         u_data <- cbind(seq(1, tmax, by=deltat), s$u);
#         colnames(u_data) <- c("time", names);
#         xmap <- ccm(u_data, E = E, lib_column = names[i], 
#                     target_column = names[j], lib_sizes = lib_sizes,
#                     num_samples = num_samples, random_libs = TRUE, replace = TRUE);
#         xmap_means <- ccm_means(xmap)
#         ccm_rho_graphs[i, j,, t] <- pmax(0, xmap_means$rho)
#         
#         plot(xmap_means$lib_size, pmax(0, xmap_means$rho), type = "l", col = "red", xlab = "Library Size", 
#              ylab = "Cross Map Skill (rho)", ylim = c(0, 1), main = paste(i, j))
#         legend(x = "topleft", legend = paste(i, "xmap", j),
#                col = c("red"), lwd = 1, bty = "n", inset = 0.02, cex = 0.8)
#       }
#       print(paste("xmap", i, j))
#     }
#   }
#   
#   return(list("ccm_rho_graphs" = ccm_rho_graphs, "svd_data" = svd_data))
# }
# 
# 
# get_adj_hankel <- function(ccm_rho_graphs, svd_data) {
#   # ccm_rho_graphs is a 4 dimensional array where the first and second
#   # dimensions represent the nodes in the system, the third
#   # dimension is the number of library sizes tried, and the fourth
#   # dimension is the number of simulation trials ran.
#   n <- dim(svd_data[[1]]$v)[1]
#   rank <- dim(ccm_rho_graphs)[1]
#   num_trials <- dim(ccm_rho_graphs)[4]
#   
#   adj_u_eig_opt_rho <- array(0, dim = c(rank, rank, num_trials))
#   max_rho <- apply(ccm_rho_graphs, c(1, 2, 4), max)
#   
#   # Search through rho values and find the one that represents 95%
#   # of the correlation maximum for each graph.
#   for (i in 1 : rank) {
#     for (j in 1 : rank) {
#       for (t in 1 : num_trials) {
#         idx <- which(ccm_rho_graphs[i, j,, t] >= 0.95 * max_rho[i, j, t])[1]
#         adj_u_eig_opt_rho[i, j, t] <- ccm_rho_graphs[i, j, idx, t]
#       }
#     }
#   }
#   
#   # Project the U eigenvalue time series adjacency matrix into the
#   # V eigenvalue basis and scale it by S where the hankel matrix H
#   # satisfies H = USV*.
#   adj <- array(0, dim = c(n, n, num_trials))
#   for (t in 1 : num_trials) {
#     s <- svd_data[[t]]
#     D <- diag(s$d[1 : rank])
#     adj[,, t] <- s$v %*% D %*% adj_u_eig_opt_rho[,, t] %*% solve(D) %*% Conj(t(s$v))
#   }
#   
#   return(apply(adj, c(1, 2), mean))
# }
