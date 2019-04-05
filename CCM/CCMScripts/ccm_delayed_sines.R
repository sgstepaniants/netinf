library(rEDM)

# Run this simulation for tmax ranging from 10 to 100. See how the rho
# curves for s1 xmap s2 and s2 xmap s1 come closer and closer to each
# other and to a correlation coefficient of 1. This implies that if we
# look closely at the initial perturbation that caused the oscillations
# of s1 (subsequently making s2 oscillate afterwords), then we have a
# better change of understanding that s1 causes s2. The farther out we
# look in time, the harder it is for us to infer the direction
# of the causal relationship. For very large tmax's, CCM returns that
# s1 and s2 cause each other equally.

# Time Vector
tmax <- 10
time <- seq(1, tmax, 0.1)
# E is optimally found using the simplex projection algorithm
E <- 4

trials <- 10
libs <- seq(10, length(time)-10, by = 10)
noise <- 0.01

# Generate Crossmap Between Two Delayed Sine Curves
a_xmap_b_rho <- numeric(length(libs))
b_xmap_a_rho <- a_xmap_b_rho
for (i in 1:trials) {
  time1 <- time
  time1[time1 < pi] <- 0
  time2 <- time
  time2[time2 < 2 * pi] <- 0
  
  # s1 begins oscillating earlier than s2
  nncoupled_model <- data.frame(time, s1=sin(time1), s2=sin(time2 - 2))
  nncoupled_model <- nncoupled_model + cbind(0,matrix(rnorm(n*tmax,mean=0,sd=noise),tmax,n))
  
  a_xmap_b <- ccm(nncoupled_model, E = E, lib_column = "s1", 
                  target_column = "s2", lib_sizes = libs,
                  num_samples = 100, random_libs = TRUE, replace = TRUE);
  b_xmap_a <- ccm(nncoupled_model, E = E, lib_column = "s2",
                  target_column = "s1", lib_sizes = libs,
                  num_samples = 100, random_libs = TRUE, replace = TRUE);
  
  a_xmap_b_means <- ccm_means(a_xmap_b);
  b_xmap_a_means <- ccm_means(b_xmap_a);
  a_xmap_b_rho <- a_xmap_b_rho + a_xmap_b_means$rho
  b_xmap_a_rho <- b_xmap_a_rho + b_xmap_a_means$rho
  print(i)
}
a_xmap_b_rho <- a_xmap_b_rho / trials
b_xmap_a_rho <- b_xmap_a_rho / trials

par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0))
ab <- pmax(0, a_xmap_b_rho);
ba <- pmax(0, b_xmap_a_rho);

plot(a_xmap_b_means$lib_size, ab, type = "l", col = "red", xlab = "Library Size", 
     ylab = "Cross Map Skill (rho)", ylim = c(0, 1))
lines(b_xmap_a_means$lib_size, ba, col = "blue")
legend(x = "topleft", legend = c("s1 xmap s2", "s2 xmap s1"),
       col = c("red", "blue"), lwd = 1, bty = "n", inset = 0.02, cex = 0.8)
