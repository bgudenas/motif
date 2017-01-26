library(GA)
library(seqLogo)

source("./src/Sequence_Functions.R")
Seqs = Motif_Implanter(seqs = Rando_DNA_Strings(seq_num = 200, width = 500), motif ="GGGGCTTAC", var_num = 1 )
background = Motif_Implanter(seqs = Rando_DNA_Strings(seq_num = 200, width = 500), motif ="", var_num = 0 )

Seqs =c(Seqs, background)
positives = 200
motif_width = 9
# 
# priors = letterFrequency(DNAStringSet(Seqs), letters = DNA_BASES)
# prior.params = colSums(priors)/sum(priors)

source("./src/PWM_functions.R")
GA = ga("binary", nBits = (4*motif_width),
      fitness = PWM_fitness, population = seeder,
      crossover = PWM_crossover, mutation = PWM_mutation,
      popSize=300, maxiter=400,
      pmutation= 0.5, pcrossover=0.6, run = 50, monitor = plot)


pwm = matrix(GA@solution, nrow=4)
rownames(pwm) = c("A", "C", "G", "T")
seqLogo(pwm)


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
