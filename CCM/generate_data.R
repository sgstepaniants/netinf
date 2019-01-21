library(deSolve)

two_species_data <- function(init_conds, num_iter, beta=0.2) {
  time <- 1:num_iter
  df <- data.frame(time, x = 0, y = 0)
  df[1,] = c(1, init_conds)
  for (t in 1:(num_iter-1)) {
    df[t+1,'x'] <- 3.9*df[t,'x']*(1-df[t,'x']-beta*df[t,'y'])
    df[t+1,'y'] <- 3.7*df[t,'y']*(1-df[t,'y']-0.2*df[t,'x'])
  }
  
  return(df)
}

five_species_data <- function(init_conds, num_iter) {
  time <- 1:num_iter
  df <- data.frame(time, y1 = 0, y2 = 0, y3 = 0, y4 = 0, y5 = 0)
  df[1,] = c(1, init_conds)
  for (t in 1:(num_iter-1)) {
    df[t+1,'y1'] <- df[t,'y1']*(4-4*df[t,'y1']-2*df[t,'y2']-0.4*df[t,'y3'])
    df[t+1,'y2'] <- df[t,'y2']*(3.1-0.31*df[t,'y1']-3.1*df[t,'y2']-0.93*df[t,'y3'])
    df[t+1,'y3'] <- df[t,'y3']*(2.12-0.636*df[t,'y1']+0.636*df[t,'y2']-2.12*df[t,'y3'])
    df[t+1,'y4'] <- df[t,'y4']*(3.8-0.111*df[t,'y1']-0.011*df[t,'y2']+0.131*df[t,'y3']-3.8*df[t,'y4'])
    df[t+1,'y5'] <- df[t,'y5']*(4.1-0.082*df[t,'y1']-0.111*df[t,'y2']+0.125*df[t,'y3']-4.1*df[t,'y5'])
  }
  
  return(df)
}


# Nearest Neighbor Spring Mass Model
# global variable to keep track of past node positions
prev_pos <- matrix()
thresh <- 0.1
past <- 5
wait <- 5
last_root <- 1 - wait

# K is the matrix of spring constants but it is also the adjacency
# matrix for the network.
nncoupled_data <- function(n, time, K, p_fn, v_fn, m_fn, c_fn, bc="circ", pert=0) {
  pos <- p_fn(n, bc)
  vel <- v_fn(n)
  
  mass <- m_fn(n)
  damp <- c_fn(n)
  
  state <- c(pos=pos,  # initial positions
             vel=vel)  # initial velocities
  
  prev_pos <<- matrix(pos, n, 1)
  last_root <<- 1
  
  bc_num <- 0
  if (bc == "circ") {
    bc_num <- 1
  } else if (bc == "fixed") {
    bc_num <- 2
  }
  
  parameters <- c(n=n,           # number of nodes
                  m=mass,        # masses
                  K=K,           # matrix of spring constants (also functions as adjacency matrix)
                  c=damp,        # damping coefficients
                  bc_num=bc_num, # boundary conditions
                  pert=pert)     # perturbation strength
  
  out = lsode(y=state, times=time, func=odeNN, parms=parameters,
              events=list(func=event, root=TRUE), rootfun=root);
  out <- out[, c("time", paste("pos", 1:n, sep=""))]
  return(out)
}

root <- function(t, y, parms) {
  with(as.list(parms), {
    num <- 1
    end <- ncol(prev_pos)
    if (end > past && t >= last_root + wait) {
      amp <- rowrange(prev_pos[, (end - past + 1):end]) / 2
      num <- 1 - any(amp < thresh)
      
      # mess with the recorded previous node positions to ensure that
      # the nodes are not kicked on two consecutive iterations
      prev_pos[, end - 1] <- prev_pos[, end] - 4 * amp
    }
    
    if (num == 0) {
      last_root <<- t
    }
    return(num)
  })
}

event <- function(t, y, parms) {
  with(as.list(parms), {
    # If some of the nodes have amplitudes smaller than a certain threshold,
    # kick them
    if (n > 1) {
      m <- parms[paste("m", 1:n, sep="")]
    }
    
    if (ncol(prev_pos) > past) {
      amp <- rowrange(prev_pos[, (ncol(prev_pos) - past + 1): ncol(prev_pos)]) / 2
      y[(n+1):(2*n)] <- y[(n+1):(2*n)] - pert * rnorm(n) * (amp < thresh) / m;
    }
    return(y)
  })
}

odeNN <- function(t, state, parms) {
  with(as.list(parms), {
    p <- state["pos"]
    v <- state["vel"]
    K <- parms["K"]
    if (n > 1) {
      p <- state[grep('^pos', names(state))]
      v <- state[grep('^vel', names(state))]
      m <- parms[grep('^m', names(parms))]
      
      K_flat <- parms[grep('^K', names(parms))]
      num <- sqrt(length(K_flat))
      K <- matrix(K_flat, nrow = num, ncol = num)
      
      c <- parms[grep('^c', names(parms))]
    }
    
    # update global variable of all previous node positions
    prev_pos <<- cbind(prev_pos, p)
    
    if (bc_num == 1) { # circular bc's
      if (nrow(K) == n && ncol(K) == n) {
        p <- p - unifpfn(n, "circ")
        r <- matrix(p, nrow=n, ncol=n)
        dist <- t(r) - r
        force <- rowSums((K + t(K)) * dist)
      } else {
        stop("spring constants for circular boundary conditions must have dimension n x n")
      }
    } else if (bc_num == 2) { # fixed bc's (wall anchored)
      if (nrow(K) == n + 2 && ncol(K) == n + 2) {
        p <- p - unifpfn(n, "fixed")
        p <- c(0, p, 0)
        r <- matrix(p, nrow=n + 2, ncol=n + 2)
        dist <- t(r) - r
        
        force <- rowSums(K * dist)
        force <- force[2 : (n + 1)]
      } else {
        stop('spring constants for fixed boundary conditions must have dimension (n + 2) x (n + 2)')
      }
    } else { # free bc's (no springs on ends)
      if (nrow(K) == n && ncol(K) == n) {
        p <- p - unifpfn(n, "free")
        r <- matrix(p, nrow=n, ncol=n)
        dist <- t(r) - r
        force <- rowSums(K * dist)
      } else {
       stop('spring constants for free boundary conditions must have dimension n x n')
     }
    }
    
    # calculate the velocity and acceleration
    vel <- as.numeric(v)
    acc <- (force - c * v) / m
    
    dd <- c(dPos=vel, dVel=acc)
    return(list(dd))
  })
}

rowrange <- function(x) {
  return(apply(x, 1, function(x) diff(range(x))))
}

circshift <- function(df, n = 1) {
  dft <- t(t(df))
  df_shifted <- rbind(tail(dft,n), head(dft,-n), deparse.level = 0)
  return(t(df_shifted))
}






# Kuramoto Oscillator Model
kuramoto_data <- function(A, time, K, randicfn, randwfn) {
  n <- dim(A)[1]
  m <- length(time)
  
  pos <- randicfn(n)
  freq <- randwfn(n)
  
  state <- c(pos=pos)  # initial positions
  
  parameters <- c(n=n,     # number of nodes
                  adj=A,   # adjacency matrix
                  w=freq,  # natural frequencies
                  K=K)     # connection strength
  
  out = ode(y=state, times=time, func=odeKur, parms=parameters);
  out <- out[, c("time", paste("pos", 1:n, sep=""))]
  return(out)
}

odeKur <- function(t, state, parameters) {
  with(as.list(parameters), {
    # every entry A[i, j] ith row and jth column tells us if
    # the jth node influences the ith node
    A <- parameters["adj"]
    w <- parameters["w"]
    if (n > 1) {
      adj <- parameters[paste("adj", 1:(n*n), sep="")]
      A <- matrix(adj, nrow = n, ncol = n)
      w <- parameters[paste("w", 1:n, sep="")]
    }
    
    y <- state
    r <- matrix(y, nrow=n, ncol=n)
    dy <- w + K/n*rowSums(A * sin(t(r)-r))
    return(list(dy))
  })
}
