library(rEDM)

# Number of Nodes
n = 2
A = matrix(c(0, 1, 0, 0), nrow = n, ncol = n)
# Time Vector
tmax = 30
time = 1:tmax
E = 10
K = 1

trials = 100

# Generate Nearest-Neighbor Spring Mass Model Data
a_xmap_b_rho <- numeric(tmax)
b_xmap_a_rho <- numeric(tmax)
for (i in 1:trials) {
  kuramoto_model <- kuramoto_data(A, time, K, randicfn, randwfn)
  kuramoto_model <- cos(kuramoto_model)
  
  a_xmap_b <- ccm(kuramoto_model, E = E, lib_column = "pos1", 
                  target_column = "pos2", lib_sizes = seq(10, 500, by = 10),
                  num_samples = 100, random_libs = TRUE, replace = TRUE);
  b_xmap_a <- ccm(kuramoto_model, E = E, lib_column = "pos2",
                  target_column = "pos1", lib_sizes = seq(10, 500, by = 10),
                  num_samples = 100, random_libs = TRUE, replace = TRUE);
  
  a_xmap_b_means <- ccm_means(a_xmap_b);
  b_xmap_a_means <- ccm_means(b_xmap_a);
  a_xmap_b_rho <- a_xmap_b_rho + a_xmap_b_means$rho
  b_xmap_a_rho <- b_xmap_a_rho + b_xmap_a_means$rho
}
a_xmap_b_rho <- a_xmap_b_rho / trials
b_xmap_a_rho <- b_xmap_a_rho / trials

par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0))
ab <- pmax(0, a_xmap_b_rho);
ba <- pmax(0, b_xmap_a_rho);

plot(a_xmap_b_means$lib_size, ab, type = "l", col = "red", xlab = "Library Size", 
     ylab = "Cross Map Skill (rho)", ylim = c(0, 1), main = "a = 1, b = 2")
lines(b_xmap_a_means$lib_size, ba, col = "blue")
legend(x = "topleft", legend = c("a xmap b", "b xmap a"),
       col = c("red", "blue"), lwd = 1, bty = "n", inset = 0.02, cex = 0.8)
