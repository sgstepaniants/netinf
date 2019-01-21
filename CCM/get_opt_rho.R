get_adj <- function(ccm_rho_graphs) {
  # ccm_rho_graphs is a 4 dimensional array where the first and second
  # dimensions represent the nodes in the system, the third
  # dimension is the number of library sizes tried, and the fourth
  # dimension is the number of simulation trials ran.
  n <- dim(ccm_rho_graphs)[1]
  trials <- dim(ccm_rho_graphs)[4]
  
  adj_opt_rho <- array(0, dim = c(n, n, trials))
  max_rho <- apply(ccm_rho_graphs, c(1, 2, 4), max)
  
  # Search through rho values and find the one that represents 95%
  # of the correlation maximum for each graph.
  for (i in 1 : n) {
    for (j in 1 : n) {
      for (t in 1 : trials) {
        idx <- which(ccm_rho_graphs[i, j,, t] >= 0.95 * max_rho[i, j, t])[1]
        adj_opt_rho[i, j, t] <- ccm_rho_graphs[i, j, idx, t]
      }
    }
  }
  
  return(adj_opt_rho)
}

get_adj_hankel <- function(ccm_rho_graphs, svd_data) {
  # ccm_rho_graphs is a 4 dimensional array where the first and second
  # dimensions represent the nodes in the system, the third
  # dimension is the number of library sizes tried, and the fourth
  # dimension is the number of simulation trials ran.
  n <- dim(svd_data[[1]]$v)[1]
  rank <- dim(ccm_rho_graphs)[1]
  trials <- dim(ccm_rho_graphs)[4]
  
  adj_u_eig_opt_rho <- array(0, dim = c(rank, rank, trials))
  max_rho <- apply(ccm_rho_graphs, c(1, 2, 4), max)
  
  # Search through rho values and find the one that represents 95%
  # of the correlation maximum for each graph.
  for (i in 1 : rank) {
    for (j in 1 : rank) {
      for (t in 1 : trials) {
        idx <- which(ccm_rho_graphs[i, j,, t] >= 0.95 * max_rho[i, j, t])[1]
        adj_u_eig_opt_rho[i, j, t] <- ccm_rho_graphs[i, j, idx, t]
      }
    }
  }
  
  # Project the U eigenvalue time series adjacency matrix into the
  # V eigenvalue basis and scale it by S where the hankel matrix H
  # satisfies H = USV*.
  adj <- array(0, dim = c(n, n, trials))
  for (t in 1 : trials) {
    s <- svd_data[[t]]
    D <- diag(s$d[1 : rank])
    adj[,, t] <- s$v %*% D %*% adj_u_eig_opt_rho[,, t] %*% solve(D) %*% Conj(t(s$v))
  }
  
  return(apply(adj, c(1, 2), mean))
}
