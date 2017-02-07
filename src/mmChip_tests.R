# Motif GA test on mouse Chip-seq data ------------------------------------------
library(stringr)
library(Biostrings)
library(biomaRt)
source("http://www.bioconductor.org/biocLite.R")
biocLite("seqLogo")
biocLite( "BSgenome.Mmusculus.UCSC.mm8")
library(BSgenome)
library(BSgenome.Mmusculus.UCSC.mm8)
l


mmChip = read.csv("./Data/mmc3.csv")

Sequences <- data.frame(lapply(mmChip, as.character), stringsAsFactors=FALSE)

for (col in 1:ncol(mmChip)){
        
    for ( i in 1:length(mmChip[ ,col])) {
        if ( mmChip[i, col] != ""){
        peak = str_split(mmChip[i, col] , "-")[[1]]
        chrom = str_split(peak[[1]][1], ":")[[1]][1]
        start = as.numeric(str_split(peak[1], ":")[[1]][2])
        end = as.numeric(peak[2])
        center = median(start:end)
        
        genome_sub = Mmusculus[[chrom]][(center-50):(center+49)]
        Sequences[i, col] = as.character(genome_sub)
        
        }
    }
}


write.csv("./Data/mmChip_Sequences.csv", x = Sequences)


# test PWM_GA on mmchip data ----------------------------------------------
Sequences = read.csv("./Data/mmChip_Sequences.csv", stringsAsFactors = FALSE)

library(GA)
library(seqLogo)

source("./src/Sequence_Functions.R")
# Seqs = Motif_Implanter(seqs = Rando_DNA_Strings(seq_num = 200, width = 500), motif ="ATGCACAA", var_num = 1 )
# background = Motif_Implanter(seqs = Rando_DNA_Strings(seq_num = 200, width = 500), motif ="", var_num = 0 )


Seqs = Sequences$c.Myc.bound.loci
Seqs = Seqs[Seqs!=""]
positives = length(Seqs)

background = Seqs
for ( i in 1:length(background)){
    background[i] = dinuc_shuffle(background[i])
}
# background = Sequences$Sox2.bound.loci
# background = background[background!=""]
# # background = background[1:length(Seqs)]

Seqs =c(Seqs, background)
motif_width = 8
# priors = letterFrequency(DNAStringSet(Seqs), letters = DNA_BASES)
# prior.params = colSums(priors)/sum(priors)

source("./src/PWM_functions.R")
GA = ga("binary", nBits = (4*motif_width),
        fitness = PWM_fitness, population = seeder,
        crossover = PWM_crossover, mutation = PWM_mutation,
        popSize= 200, maxiter=300,
        pmutation= 0.4, pcrossover=0.5, run = 20)


pwm = matrix(GA@solution, nrow=4)
rownames(pwm) = c("A", "C", "G", "T")
seqLogo(pwm)


pop = GA@population[order(GA@fitness, decreasing = TRUE), ]

pwm2 = matrix(pop[5,] , nrow=4)
rownames(pwm2) = c("A", "C", "G", "T")
seqLogo(pwm2)



vec = vector(length = length(Seqs))
for (i in 1:length(Seqs)){
    vec[i] = countPattern(pattern = "CACGTG", subject = Seqs[i])
}
