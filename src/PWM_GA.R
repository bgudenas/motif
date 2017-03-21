library(GA)
library(seqLogo)

source("./src/Sequence_Functions.R")
Seqs = Motif_Implanter(seqs = Rando_DNA_Strings(seq_num = 100, width = 200), motif ="GGGGCTTAAG", var_num = 2 )
background = Motif_Implanter(seqs = Rando_DNA_Strings(seq_num = 100, width = 200), motif ="", var_num = 0 )

Seqs =c(Seqs, background)
positives = 100
motif_width = 10
# 
# priors = letterFrequency(DNAStringSet(Seqs), letters = DNA_BASES)
# prior.params = colSums(priors)/sum(priors)

source("./src/PWM_functions.R")
GA = ga("binary", nBits = (4*motif_width),
      fitness = PWM_fitness, population = seeder,
      crossover = PWM_crossover, mutation = PWM_mutation,
      popSize=200, maxiter=800,
      pmutation= 0.3, pcrossover=0.5, run = 20, monitor = plot)


pwm = matrix(GA@solution, nrow=4)
rownames(pwm) = c("A", "C", "G", "T")
seqLogo(pwm)



PWMscoreStartingAt(pwm, org_motif)/maxScore(pwm) * 100
