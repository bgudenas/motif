
# Synthetic motif_GA tests ----------------------------------------
res_mat = matrix(data=0, nrow = 5, ncol = 8)
rownames(res_mat) = c(100, 200, 400, 800, 1600)

library(GA)
library(seqLogo)

source("./src/Sequence_Functions.R")
org_motif = "ATGGTGCTCGCATTA"
Seqs = Motif_Implanter(seqs = Rando_DNA_Strings(seq_num = 100, width = 1600), motif = "ATGGTGCTCGCATTA", var_num = 4 )
background = Motif_Implanter(seqs = Rando_DNA_Strings(seq_num = 100, width = 1600), motif ="", var_num = 0 )
positives = length(Seqs)

# Seqs = Sequences$Zfx.bound.loci
# Seqs = Seqs[Seqs!=""]
# positives = length(Seqs)
# 
# background = Seqs
# for ( i in 1:length(background)){
#     background[i] = dinuc_shuffle(background[i])
# }
# # background = Sequences$Sox2.bound.loci
# # background = background[background!=""]
# # # background = background[1:length(Seqs)]

Seqs =c(Seqs, background)
motif_width = 15
# priors = letterFrequency(DNAStringSet(Seqs), letters = DNA_BASES)
# prior.params = colSums(priors)/sum(priors)
answers = vector(mode = "numeric", length = 8)
for ( a in 1:length(answers)){

source("./src/PWM_functions.R")
GA = ga("binary", nBits = (4*motif_width),
        fitness = PWM_fitness, population = seeder,
        crossover = PWM_crossover, mutation = PWM_mutation,
        popSize= 100, maxiter=300,
        pmutation= 0.4, pcrossover=0.6, run = 20)


pwm = matrix(GA@solution, nrow=4)
rownames(pwm) = c("A", "C", "G", "T")
seqLogo(pwm)



answers[a] = 1 - (stringdist(org_motif, pwm_2_seq(pwm)) / 15)
}



res_mat[4, ] = answers


