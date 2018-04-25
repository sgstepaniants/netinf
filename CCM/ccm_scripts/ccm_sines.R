library(rEDM)

# Time Vector
tmax <- 100
time <- 1:tmax

trials <- 10
libs <- seq(10, tmax-10, by = 10)

# Generate Crossmap Between Two Sine Curves
a_xmap_b_rho <- numeric(length(libs))
b_xmap_a_rho <- a_xmap_b_rho
for (i in 1:trials) {
  # the sine that oscillates faster seems to "CCM cause" the slower sine
  nncoupled_model <- data.frame(time, s1=sin(2*time + 12), s2=50*sin(time))
  nncoupled_model <- nncoupled_model + cbind(0,matrix(rnorm(n*tmax,mean=0,sd=0.01),tmax,n))
  
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