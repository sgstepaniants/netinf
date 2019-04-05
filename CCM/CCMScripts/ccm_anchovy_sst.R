library(rEDM)

data(sardine_anchovy_sst)
anchovy_xmap_sst <- ccm(sardine_anchovy_sst, E = 3, lib_column = "anchovy", 
                        target_column = "np_sst", lib_sizes = seq(10, 80, by = 10), num_samples = 100, 
                        random_libs = TRUE, replace = TRUE)
sst_xmap_anchovy <- ccm(sardine_anchovy_sst, E = 3, lib_column = "np_sst", target_column = "anchovy", 
                        lib_sizes = seq(10, 80, by = 10), num_samples = 100, random_libs = TRUE, 
                        replace = TRUE)

a_xmap_t_means <- ccm_means(anchovy_xmap_sst)
t_xmap_a_means <- ccm_means(sst_xmap_anchovy)

par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0))
y1 <- pmax(0, a_xmap_t_means$rho)
y2 <- pmax(0, t_xmap_a_means$rho)

plot(a_xmap_t_means$lib_size, y1, type = "l", col = "red", xlab = "Library Size", 
     ylab = "Cross Map Skill (rho)", ylim = c(0, 0.25))
lines(t_xmap_a_means$lib_size, y2, col = "blue")
legend(x = "topleft", legend = c("anchovy xmap SST", "SST xmap anchovy"),
       col = c("red", "blue"), lwd = 1, bty = "n", inset = 0.02, cex = 0.8)