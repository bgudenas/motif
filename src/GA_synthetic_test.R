library(GA)
library(seqLogo)

source("./src/Sequence_Functions.R")
org_motif = "CCCGGGTAAAGTCCC"
org_motif
# Synthetic motif_GA tests ----------------------------------------
res_mat = matrix(data=0, nrow = 5, ncol = 11)
rownames(res_mat) = c( 200, 400, 600, 800, 1000)   #400, 800, 1600)


for ( i in as.numeric(rownames(res_mat))){
Seqs = Motif_Implanter(seqs = Rando_DNA_Strings(seq_num = 100, width = i), motif = org_motif, var_num = 3 )
background = Motif_Implanter(seqs = Rando_DNA_Strings(seq_num = 100, width = i), motif ="", var_num = 0 )
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
            popSize= 500, maxiter=1000,
            pmutation= 0.3, pcrossover=0.6, run = 40)
    
    
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
res_mat = t(res_mat)
rownames(res_mat) = 1:nrow(res_mat)
res_mat %>% 
    group_by(as.factor(colnames(res_mat))) %>% 
    ggplot( aes(x=rownames(res_mat), y= colnames(res_mat)), fill= colnames(res_mat)) +
    geom_boxplot() +
    facet_grid(.~colnames(res_mat))

write.csv(res_mat, "./Data/synth_results.csv")

res_mat[2, ] = answers

library(ggplot2)
library(dplyr)

means = rowMeans(res_mat)
means = c("100"=76.40431, means)

par(mar=c(5.1,5.1,4.1,2.1))
barplot(means, ylim =range(0,100), xlab="Sequence Length", ylab = "Mean Similarity", names.arg = names(means), cex.axis = 2, cex.names = 2.3, cex.lab = 2.3)
title("Motif Finding Accuracy", cex.main=2.2)
abline(h = 60, col="red",lwd =3)
abline(h = 80,lwd =3, lty=2)

text(x=5,y=82, "Noise added per Motif", cex =2)
text(x=means ,y=20, c(".15",".075","0"), cex =1)


library(ggplot2)
means =rowMeans(res_mat)
means[,1] = as.numeric(means[,1])
vals = cbind(means, names(means))
ggplot( aes(x=means, y = names(means))) +
    geom_bar(stat = "identity")

ggplot(data = as.data.frame(res_mat), aes(x=100, y= ))+
    + geom_boxplot()

par(mfrow=c(1,5))
boxplot(res_mat[,1], ylim=range(0,100))
abline(h=60, col="red",lwd=2)
boxplot(res_mat[,2], ylim=range(0,100))
abline(h=60, col="red",lwd=2)
boxplot(res_mat[,3], ylim=range(0,100))
abline(h=60, col="red",lwd=2)
boxplot(res_mat[,4], ylim=range(0,100))
abline(h=60, col="red",lwd=2)
boxplot(res_mat[,5], ylim=range(0,100))
abline(h=60, col="red",lwd=2)
