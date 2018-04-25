›4library(plotrix)

# Number of Nodes
n = 5;
# Time Vector
time = 1:500
# Boundary Conditions
bc = "circ"

singlepfn <- function(n) {
  return(c(0.1, seq(from=1/n, to=(n-1)/n, by=1/n))) # circ
  #return(c(0, seq(from=2/(n+1), to=n/(n+1), by=1/(n+1)))) # fixed
}
unifpfn <- function(n) {
  return(seq(from=0, to=(n-1)/n, by=1/n)) # circ
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

nncoupled_model <- nncoupled_data(n, time, singlepfn, zerovfn, constmfn, constkfn, bc=bc)

for (i in 1:length(time)) {
  t <- time[i]
  if (bc=="circ") {
    plot(-2:2,-2:2,type="n",main="Spring Mass Simulation")
    for (x in 1:n) {
        pos <- 2*pi*nncoupled_model[t, x+1];
        r <- 0.1;
        draw.circle(cos(pos), sin(pos), r, col=heat.colors(n)[x])
    }
  } else {
    plot(-0.25:1.25,-0.75:0.75,type="n",main="Spring Mass Simulation")
    for (x in 1:n) {
      pos <- nncoupled_model[t, x+1];
      r <- 0.1;
      draw.circle(pos, 0, r, col=heat.colors(n)[x])
    }
  }
  Sys.sleep(0.05)
  
  if (i < length(time)) {
    #dev.off()
    #plot.new()
  }
}
