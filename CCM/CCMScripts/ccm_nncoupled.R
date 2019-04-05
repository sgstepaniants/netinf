library(rEDM)

# Number of Nodes
n <- 20
E <- 41
# Time Vector
time <- 1:500
# Boundary Conditions
bc <- "circ"

# Generate Nearest-Neighbor Spring Mass Model Data
nncoupled_model <- nncoupled_data(n, time, smallrandpfn, zerovfn, constmfn, constkfn, bc=bc)

a_xmap_b <- ccm(nncoupled_model, E = E, lib_column = "pos1", 
                target_column = "pos2", lib_sizes = seq(10, 500, by = 10),
                num_samples = 100, random_libs = TRUE, replace = TRUE);
b_xmap_a <- ccm(nncoupled_model, E = E, lib_column = "pos2",
                target_column = "pos1", lib_sizes = seq(10, 500, by = 10),
                num_samples = 100, random_libs = TRUE, replace = TRUE);

a_xmap_b_means <- ccm_means(a_xmap_b);
b_xmap_a_means <- ccm_means(b_xmap_a);

par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0))
ab <- pmax(0, a_xmap_b_means$rho);
ba <- pmax(0, b_xmap_a_means$rho);

plot(a_xmap_b_means$lib_size, ab, type = "l", col = "red", xlab = "Library Size", 
     ylab = "Cross Map Skill (rho)", ylim = c(0, 1), main = "a = 1, b = 6, n = 20, E = 41")
lines(b_xmap_a_means$lib_size, ba, col = "blue")
legend(x = "topleft", legend = c("a xmap b", "b xmap a"),
       col = c("red", "blue"), lwd = 1, bty = "n", inset = 0.02, cex = 0.8)
