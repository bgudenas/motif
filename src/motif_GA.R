library(GA)

source("./src/Sequence_Functions.R")
test_seqs = Motif_Implanter(seqs = Rando_DNA_Strings(seq_num = 20, width = 30), motif ="GGGGGGGG" )
background = Motif_Implanter(seqs = Rando_DNA_Strings(seq_num = 20, width = 30), motif ="" )

Seqs = c(test_seqs, background)
names(Seqs)[1:20] = "Positive"
names(Seqs)[21:40] = "Negative"

pos = 20
motif_width = 8

pos.freq = letterFrequency(DNAStringSet(Seqs), letters = DNA_BASES)[1:pos, ]
pos.params = colSums(pos.freq)/sum(pos.freq)

neg.freq = letterFrequency(DNAStringSet(Seqs), letters = DNA_BASES)[(pos+1):length(Seqs), ]
neg.params = colSums(neg.freq)/sum(neg.freq)

GA=ga("real-valued", min = rep(1, length(Seqs)), max = unlist(lapply(Seqs, function(x) nchar(x)-motif_width)),
      fitness = PWM_dist, population = seeder,
      # fitness=itness, crossover = FB_crossover,population = seeder,
      popSize=1000, maxiter=1,
      pmutation= 0.3, pcrossover=0.7,run=1)










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
