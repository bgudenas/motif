library(GA)
library(seqLogo)
source("./src/Sequence_Functions.R")

gencode = read.csv("./Data/gencode_seqs.csv")
colnames(gencode)[ncol(gencode)] = "Seqs"

Seqs = as.character(gencode$Seqs[gencode$RPKM1.RPKM2 > 1])
background = as.character(gencode$Seqs[gencode$RPKM1.RPKM2 < 1 ])
positives = length(Seqs)
motif_width = 14

Seqs = c(Seqs, background)
  
pwm = matrix(GA@solution, nrow=4)
rownames(pwm) = c("A", "C", "G", "T")
seqLogo(pwm)
pwm = pwm + 0.0001

vec =vector("numeric", length = length(Seqs))
for ( i in 1:length(Seqs)){
    vec[i] = countPWM(pwm, Seqs[i], min.score = "60%")
}


sum(vec[1:positives]) / sum(vec[(positives+1):length(vec)])
sum(vec[(positives+1):length(vec)])


lncrna_lens = vector(mode="numeric",length = length(Seqs))
for ( i in 1:length(Seqs)){
    lncrna_lens[i] = nchar(Seqs[i])
    
}


sum(lncrna_lens)


background = Motif_Implanter(seqs = Rando_DNA_Strings(seq_num = 1749, width = 1400), motif ="", var_num = 0 )
bvec =vector("numeric", length = length(background))
for ( i in 1:length(background)){
    bvec[i] = countPWM(pwm, background[i], min.score = "60%")
}
sum(bvec[1:positives]) / sum(vec[(positives+1):length(vec)])
