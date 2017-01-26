library(GA)
library(seqLogo)

source("./src/Sequence_Functions.R")
Seqs = Motif_Implanter(seqs = Rando_DNA_Strings(seq_num = 400, width = 300), motif ="GGGGGGGG" )
background = Motif_Implanter(seqs = Rando_DNA_Strings(seq_num = 200, width = 300), motif ="" )
# 
# Seqs = c(test_seqs, background)
# names(Seqs)[1:200] = "Positive"
# names(Seqs)[201:400] = "Negative"

positives = 400
motif_width = 8

priors = letterFrequency(DNAStringSet(Seqs), letters = DNA_BASES)
prior.params = colSums(priors)/sum(priors)


GA=ga("binary", nBits = length(Seqs),
      fitness = PWM_dist, population = seeder,
      crossover = pt_crossover, mutation = mutation,
      popSize=400, maxiter=200,
      pmutation= 0.5, pcrossover=0.6, run = 25, monitor = plot)


motif = DNAStringSet(Seqs, start = GA@solution, end = GA@solution + (motif_width-1))
cons = consensusMatrix(motif, as.prob = TRUE)[1:4, ]


pwm = PWM(motif)
seqLogo(cons)


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
