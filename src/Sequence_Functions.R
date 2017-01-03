source("https://bioconductor.org/biocLite.R")
biocLite()
biocLite("BiocInstaller")
biocLite("BiocUpgrade")
biocLite("S4Vectors")
biocLite("Biostrings")
# Create function to generate random sequences  ---------------------------
library(Biostrings)


Rando_DNA_Strings = function(seq_num, width){
#generates a list of random DNA sequence 
    seqs = vector(mode="character", length = seq_num)
    for ( i in 1:seq_num){
    seqs[i] = paste(sample(DNA_ALPHABET[1:4], size=width, replace=TRUE), collapse = "")
    }
    return(seqs)
}


Motif_Implanter = function(seqs, motif){
    library(stringr)
    seq_length = nchar(seqs[1])
     # motif = paste(sample(DNA_ALPHABET[1:4], size = motif_length, replace=TRUE), collapse = "")
     motif_length=nchar(motif)
    for ( i in 1:length(seqs)){
        DNA = seqs[i]
        index = sample(1:(seq_length-motif_length), 1)
        prefix = str_sub(DNA, start = 1, end = index)
        suffix = str_sub(DNA, start = index + motif_length + 1, end = nchar(DNA))
        new_seq = str_c(prefix,motif,suffix)
        seqs[i] = new_seq
    }
return(seqs)
}


test_seqs = Motif_Implanter(seqs = Rando_DNA_Strings(seq_num = 20, width = 30), motif ="GGGGGGGG" )
background = Motif_Implanter(seqs = Rando_DNA_Strings(seq_num = 20, width = 30), motif ="" )

        
        


    