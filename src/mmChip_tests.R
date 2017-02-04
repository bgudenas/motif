# Motif GA test on mouse Chip-seq data ------------------------------------------
library(stringr)
library(Biostrings)
library(biomaRt)
source("http://www.bioconductor.org/biocLite.R")
biocLite( "BSgenome.Mmusculus.UCSC.mm8")
library(BSgenome)
library(BSgenome.Mmusculus.UCSC.mm8)


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



save.image("./Data/Oct4_CtCF_mm8_seqs.RData")


# test PWM_GA on mmchip data ----------------------------------------------

load(file="./Data/Oct4_CtCF_mm8_seqs.RData")
library(GA)
library(seqLogo)

source("./src/Sequence_Functions.R")
# Seqs = Motif_Implanter(seqs = Rando_DNA_Strings(seq_num = 200, width = 500), motif ="GGGGCTTAC", var_num = 1 )
# background = Motif_Implanter(seqs = Rando_DNA_Strings(seq_num = 200, width = 500), motif ="", var_num = 0 )


Seqs = Sequences$Seq
negative = CTCFSequences$Seq[sample(1:nrow(CTCFSequences), 1000)]
positives = length(Seqs)

Seqs =c(Seqs, negative)
motif_width = 8
# 
# priors = letterFrequency(DNAStringSet(Seqs), letters = DNA_BASES)
# prior.params = colSums(priors)/sum(priors)

source("./src/PWM_functions.R")
GA = ga("binary", nBits = (4*motif_width),
        fitness = PWM_fitness, population = seeder,
        crossover = PWM_crossover, mutation = PWM_mutation,
        popSize=100, maxiter=400,
        pmutation= 0.5, pcrossover=0.6, run = 50, monitor = plot)


pwm = matrix(GA@solution, nrow=4)
rownames(pwm) = c("A", "C", "G", "T")
seqLogo(pwm)
