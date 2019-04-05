library(rEDM)

#data(two_species_model)
two_species_model <- two_species_data(c(0.2, 0.2), 1000)
x_xmap_y <- ccm(two_species_model, E = 3, lib_column = "x", 
                        target_column = "y", lib_sizes = seq(10, 80, by = 10),
                        num_samples = 100, random_libs = TRUE, replace = TRUE);
y_xmap_x <- ccm(two_species_model, E = 3, lib_column = "y",
                        target_column = "x", lib_sizes = seq(10, 80, by = 10),
                        num_samples = 100, random_libs = TRUE, replace = TRUE);

x_xmap_y_means <- ccm_means(x_xmap_y);
y_xmap_x_means <- ccm_means(y_xmap_x);

par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0))
y1 <- pmax(0, x_xmap_y_means$rho);
y2 <- pmax(0, y_xmap_x_means$rho);

plot(x_xmap_y_means$lib_size, y1, type = "l", col = "red", xlab = "Library Size", 
     ylab = "Cross Map Skill (rho)", ylim = c(0, 1))
lines(y_xmap_x_means$lib_size, y2, col = "blue")
legend(x = "topleft", legend = c("x xmap y", "y xmap x"),
       col = c("red", "blue"), lwd = 1, bty = "n", inset = 0.02, cex = 0.8)
