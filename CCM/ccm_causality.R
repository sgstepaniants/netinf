library(rEDM)
library(lattice)

rgb.palette <- colorRampPalette(c("blue", "red"), space = "rgb")
#par(mfrow = c(2, 2))

# ANCHOVY SST DATA
data(sardine_anchovy_sst)
names <- c('anchovy', 'np_sst')
adj = get_ccm_means(sardine_anchovy_sst, 3, names)
print(adj)

# correlation coefficient threshold
thresh <- 0.1
mat <- apply(adj, c(1,2), max)
mat[mat < thresh] <- 0
mat[mat >= thresh] <- 1

levelplot(t(mat), main="Adjacency Matrix (anchovy, sst)",
          xlab="", ylab="", ylim=c(length(names) + 0.5, 0.5),
          col.regions=rgb.palette(120))

# TWO SPECIES MODEL
two_species <- two_species_data(c(0.1, 0.1, 0.1, 0.1, 0.1), 100)
names <- c('x', 'y')
adj = get_ccm_means(two_species, 3, names, lib_sizes=seq(10, 200, by = 10))
print(adj)

thresh <- 0.7
mat <- apply(adj, c(1,2), max)
mat[mat < thresh] <- 0
mat[mat >= thresh] <- 1

levelplot(t(mat), main="Adjacency Matrix Two Species",
          xlab="", ylab="", ylim=c(length(names) + 0.5, 0.5),
          col.regions=rgb.palette(120))

# FIVE SPECIES MODEL
five_species <- five_species_data(c(0.1, 0.1, 0.1, 0.1, 0.1), 100)
names <- c('y1', 'y2', 'y3', 'y4', 'y5')
adj = get_ccm_means(five_species, 3, names, lib_sizes=seq(10, 200, by = 10))

thresh <- 0.5
mat <- apply(adj, c(1,2), max)
mat[mat < thresh] <- 0
mat[mat >= thresh] <- 1

levelplot(t(mat), main="Adjacency Matrix Five Species",
          xlab="", ylab="", ylim=c(length(names) + 0.5, 0.5),
          col.regions=rgb.palette(120))






# NEAREST NEIGHBOR SPRING MODEL
# Number of Nodes
n = 3;
# Time Vector
time = 1:500
# Boundary Conditions
bc = "free"

# Initital Condition Functions
randpfn <- function(n) {
  return(runif(n, 0, 1))
}
unifpfn <- function(n) {
  return(0:1/n:(n-1)/n)
}
zerovfn <- function(n) {
  return(rep(0, n))
}
constmfn <- function(n) {
  return(rep(1, n))
}
constkfn <- function(n) {
  return(rep(1, n))
}

nncoupled_model <- nncoupled_data(n, time, randpfn, zerovfn, constmfn, constkfn, bc=bc)
names <- paste("pos", 1:n, sep="")
adj = get_ccm_means(nncoupled_model, 3, names, lib_sizes=seq(10, 200, by = 10))

thresh <- 0.5
mat <- apply(adj, c(1,2), max)
#mat[mat < thresh] <- 0
#mat[mat >= thresh] <- 1

levelplot(t(mat), main="Adjacency Matrix NNCoupled",
          xlab="", ylab="", ylim=c(length(names) + 0.5, 0.5),
          col.regions=rgb.palette(120))

