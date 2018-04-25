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

nncoupled_data <- function(n, time, p_fn, v_fn, m_fn, k_fn, c_fn, bc="circ", pert=0) {
  pos <- p_fn(n, bc)
  vel <- v_fn(n)
  
  mass <- m_fn(n)
  spring <- k_fn(n, bc)
  damp <- c_fn(n)
  
  if (length(spring) < n - 1 || length(spring) > n + 1) {
    stop('spring constants must have length n - 1, n, or n + 1')
  }
  
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
                  k=spring,      # spring constants
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

odeNN <- function(t, state, parameters) {
  with(as.list(parameters), {
    p <- state["pos"]
    v <- state["vel"]
    if (n > 1) {
      p <- state[paste("pos", 1:n, sep="")]
      v <- state[paste("vel", 1:n, sep="")]
      m <- parameters[paste("m", 1:n, sep="")]
      k <- parameters[paste("k", 1:n, sep="")]
      c <- parameters[paste("c", 1:n, sep="")]
    }
    
    p0 = circshift(p, 1)
    p1 = p
    p2 = circshift(p, -1)
    
    shift0 <- p1 - p0
    shift1 <- p2 - p1
    k0 <- c()
    k1 <- c()
    
    if (bc_num == 1) {  # n springs for circular bc's
      if (length(k) == n) {
        k0 = circshift(k, 1)
        k1 = k
      } else {
        stop("spring constants for circular boundary conditions must have length n")
      }
    
      shift0[1] = shift0[1] + 1;
      shift1[n] = shift1[n] + 1;
    } else if (bc_num == 2) {  # n + 1 springs for fixed bc's (wall anchored)
      if (length(k) == n) {
        # if n spring constants given, make the spring constants of the
        # springs attached to the walls the same
        k0 = circshift(k, 1);
        k1 = k;
      } else if (length(k) == n + 1) {
        # if n + 1 spring constants given, set the prev spring constants
        # (k0) to the first n of these and set the next spring
        # constants (k1) to the last n of these.
        k0 = head(k, n);
        k1 = tail(k, n);
      } else {
        stop('spring constants for fixed boundary conditions must have length n or n + 1')
      }
      shift0[1] = p[1];
      shift1[n] = 1 - p[n];
    } else {  # n - 1 springs for free bc's (no springs on ends)
      if (length(k) == n - 1) {
        # if n - 1 spring constants given, make outer spring constants
        # 0 (no springs)
        k0 = append(0, k);
        k1 = append(k, 0);
      } else if (length(k) == n) {
        # if n spring constants given, 0 out the last spring constant
        # because we only need n - 1 constants for free boundaries
        k[n] = 0;
        k0 = circshift(k, 1);
        k1 = k;
     } else {
       stop('spring constants for free boundary conditions must have length n - 1 or n')
     }
     shift0[1] = 0;
     shift1[n] = 0;
    }
    
    # calculate y''
    vel_thresh <- 0.1;
    acc_thresh <- 0.1;
    vel <- as.numeric(v)
    acc <- (-k0*shift0+k1*shift1-c*v)/m
    
    # update global variable of all previous node positions
    prev_pos <<- cbind(prev_pos, p)

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




# Initialization functions for spring mass nearest neighbor model
randpfn <- function(n, bc) {
  p <- runif(n, 0, 1)
  return(p)
}

smallrandpfn <- function(n, bc) {
  p <- seq(from=0, to=(n-1)/n, by=1/n) + rnorm(n, mean=0, sd=1/(2*n)) # circ
  if (bc == "fixed") {
    p <- seq(from=1/(n+1), to=n/(n+1), by=1/(n+1)) + rnorm(n, mean=0, sd=1/(2*(n+1)))
  } else if (bc == "free") {
    p <- seq(from=0, to=1, by=1/(n-1)) + rnorm(n, mean=0, sd=1/(2*(n-1)))
  }
  return(p)
}

singlepfn <- function(n, bc) {
  p <- c(-1/n, seq(from=1/n, to=(n-1)/n, by=1/n))  # circ
  if (bc == "fixed") {
    p <- c(0, seq(from=2/(n+1), to=n/(n+1), by=1/(n+1)))
  } else if (bc == "free") {
    p <- c(-1/(n-1), seq(from=1/(n-1), to=1, by=1/(n-1)))
  }
  return(p)
}

unifpfn <- function(n, bc) {
  p <- seq(from=0, to=(n-1)/n, by=1/n)  # circ
  if (bc == "fixed") {
    p <- seq(from=1/(n+1), to=n/(n+1), by=1/(n+1))
  } else if (bc == "free") {
    p <- seq(from=0, to=1, by=1/(n-1))
  }
  return(p)
}

zerovfn <- function(n) {
  return(rep(0, n))
}

constmfn <- function(n) {
  return(rep(1, n))
}

randmfn <- function(n) {
  return(runif(n, 0, 1))
}

largerandmfn <- function(n) {
  return(runif(n, 0, 100))
}

constkfn <- function(n, bc) {
  k <- rep(1, n)  # circ
  if (bc == "fixed") {
    k <- rep(1, n+1)
  } else if (bc == "free") {
    k <- rep(1, n-1)
  }
  return(k)
}

randkfn <- function(n, bc) {
  k <- runif(n, 0, 1) # circ
  if (bc == "fixed") {
    k <- runif(n+1, 0, 1)
  } else if (bc == "free") {
    k <- runif(n-1, 0, 1)
  }
  return(k)
}

largerandkfn <- function(n, bc) {
  k <- runif(n, 0, 100) # circ
  if (bc == "fixed") {
    k <- runif(n+1, 0, 100)
  } else if (bc == "free") {
    k <- runif(n-1, 0, 100)
  }
  return(k)
}

constcfn <- function(n) {
  return(rep(0.3, n))
}

zerocfn <- function(n) {
  return(rep(0, n))
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
    A <- parameters["adj"]
    w <- parameters["w"]
    if (n > 1) {
      adj <- parameters[paste("adj", 1:(n*n), sep="")]
      A <- matrix(adj, nrow = n, ncol = n)
      w <- parameters[paste("w", 1:n, sep="")]
    }
    
    y <- state
    r <- matrix(y, nrow=n, ncol=length(y))
    dy <- w + K/n*rowSums(A * sin(t(r)-r))
    return(list(dy))
  })
}

# Initialization function for Kuramoto model
randicfn <- function(n) {
  ic <- 2*pi*runif(n, 0, 1);  # uniform [0, 2pi]
  return(ic)
}

randwfn <- function(n) {
  w <- 2*runif(n, 0, 1)  # uniform [0, 2]
  return(w)
}

