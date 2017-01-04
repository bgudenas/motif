library(GA)


GA=ga("real-valued", min = rep(0, length(test_seqs)), max = unlist(lapply(test_seqs, function(x) nchar(x))),
      fitness=itness, crossover = FB_crossover, 
      mutation = mutation, popSize=1000, maxiter=100,
      pmutation= 0.3, pcrossover=0.7,run=20)


f <- function(x) abs(x)+cos(x)
curve(f, -20, 20)
fitness <- function(x) -f(x)
GA <- ga(type = "real-valued", fitness = fitness, min = -20, max = 20, popSize = 100, run = 10)
GA@population



Rastrigin <- function(x1, x2)
{20 + x1^2 + x2^2 - 10*(cos(2*pi*x1) + cos(2*pi*x2))
}
x1 <- x2 <- seq(-5.12, 5.12, by = 0.1)
f <- outer(x1, x2, Rastrigin)
persp3D(x1, x2, f, theta = 50, phi = 20)
filled.contour(x1, x2, f, color.palette = jet.colors)
