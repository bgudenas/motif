library(GA)
library(seqLogo)

source("./src/Sequence_Functions.R")
org_motif = "AATTTCTTTTCCTGT"
org_motif
# Synthetic motif_GA tests ----------------------------------------
res_mat = matrix(data=0, nrow = 5, ncol = 11)
rownames(res_mat) = c(100, 200, 400, 800, 1600)   #400, 800, 1600)


for ( i in as.numeric(rownames(res_mat))){
Seqs = Motif_Implanter(seqs = Rando_DNA_Strings(seq_num = 200, width = i), motif = org_motif, var_num = 4 )
background = Motif_Implanter(seqs = Rando_DNA_Strings(seq_num = 200, width = i), motif ="", var_num = 0 )
positives = length(Seqs)



Seqs =c(Seqs, background)
motif_width = 15
# priors = letterFrequency(DNAStringSet(Seqs), letters = DNA_BASES)
# prior.params = colSums(priors)/sum(priors)

        for ( a in 1:ncol(res_mat)){
    
    source("./src/PWM_functions.R")
    GA = ga("binary", nBits = (4*motif_width),
            fitness = PWM_fitness, population = seeder,
            crossover = PWM_crossover, mutation = PWM_mutation,
            popSize= 300, maxiter=300,
            pmutation= 0.3, pcrossover=0.5, run = 40)
    
    
    pwm = matrix(GA@solution, nrow=4)
    rownames(pwm) = c("A", "C", "G", "T")
    seqLogo(pwm)
    
    res_mat[rownames(res_mat) == as.character(i), a] = PWMscoreStartingAt(pwm, org_motif)/maxScore(pwm) * 100
        }
}

# 
# 
# score = 1
# for ( i in 1:ncol(pwm)){
#     let = str_sub(string = org_motif, start = i, end = i)
#     score = score * pwm[rownames(pwm)== let  ,i]
# }
# score/maxScore(pwm)
write.csv(res_mat, "./Data/synth_results.csv")

res_mat[2, ] = answers

library(ggplot2)
library(dplyr)
medians = vector()
> for ( i in 1:nrow(res_mat)){medians[i] = median(res_mat[i,])}

barplot(medians, ylim =range(0,100), xlab="Sequence Length (nt)", ylab = "Median Accuracy (%)", names.arg = rownames(res_mat) )
title("Motif Finding Accuracy")
text(x = rowMeans(res_mat),    paste(round(rowMeans(res_mat), 3))  , cex = 1.0)

means =rowMeans(res_mat)
means[,1] = as.numeric(means[,1])
vals = cbind(means, names(means))
ggplot( aes(x=means, y = names(rowMeans(res_mat)))) +
    geom_bar(stat = "identity")

