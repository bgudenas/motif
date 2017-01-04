
# Create function to generate random sequences  ---------------------------
library(Biostrings)
library(stringr)
library(flexclust)


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

# 
# test_seqs = Motif_Implanter(seqs = Rando_DNA_Strings(seq_num = 20, width = 30), motif ="GGGGGGGG" )
# background = Motif_Implanter(seqs = Rando_DNA_Strings(seq_num = 20, width = 30), motif ="" )


# Genetic Algorithm sub-functions -----------------------------------------
seeder = function(object) {
    pop_size = object@popSize
    max_sizes = unlist(lapply(Seqs, function(x) nchar(x)))
    seeds = matrix( ncol = length(Seqs), nrow = pop_size)
    for ( i in 1:pop_size){
        seeds[i,] = unlist(lapply(X = 1:length(Seqs), function(x) sample(1:(max_sizes[x] - motif_width), 1)))
    }
    return(seeds)
}



#TODO 
# mutation = function(object, parent){
#     solution <- as.vector(object@population[parent, ])


PWM_dist = function(individual){
    pos_set = individual[1:positives]
    neg_set = individual[(positives+1) : length(individual)]
    
    motif_pos = DNAStringSet(Seqs[1:positives], start = pos_set, end = pos_set + (motif_width-1))
    motif_neg = DNAStringSet( Seqs[(positives+1) : length(individual) ], start = neg_set, end = neg_set + (motif_width-1))
    
    PWM_pos = PWM(motif_pos, prior.params = pos.params, type="log2probratio")
    PWM_neg = PWM(motif_neg, prior.params = neg.params, type="log2probratio")
    
    fitness = sum(colSums(abs(PWM_pos - PWM_neg)))
    return(fitness)
}

    
mutation = function(object, parent){
    
     individual <- as.vector(object@population[parent, ])
     
     pick = sample(1:length(Seqs),1)
     if (pick <= positives){
         pos_set = solution[1:positives]
         pos_set = pos_set[-pick]
         
         pos_picks = 1:positives
         pos_picks = pos_picks[-pick]
         
         motif = DNAStringSet(Seqs[pos_picks], start = pos_set, end = pos_set + (motif_width-1))
         
         PWM_mat = PWM(motif, prior.params = pos.params, type="log2probratio")
         
         probs = PWMscoreStartingAt(PWM_mat, Seqs[pick], starting.at = 1:(nchar(Seqs[pick])-motif_width))
         probs = probs /sum(probs)
         
         new_pick = sample(1:(nchar(Seqs[pick])-motif_width), prob = probs, 1)
         individual[pick] = new_pick
         
     } else {
         neg_set = solution[(positives+1):length(solution)]
         neg_set = neg_set[-(pick-positives)]
         
         neg_picks = (positives+1):length(solution)
         neg_picks = neg_picks[-(pick-positives)]
         
         motif = DNAStringSet(Seqs[neg_picks], start = neg_set, end = neg_set + (motif_width-1))
         
         PWM_mat = PWM(motif, prior.params = neg.params, type="log2probratio")
         
         probs = PWMscoreStartingAt(PWM_mat, Seqs[pick], starting.at = 1:(nchar(Seqs[pick])-motif_width))
         probs = probs /sum(probs)
         
         new_pick = sample(1:(nchar(Seqs[pick])-motif_width), prob = probs, 1)
         individual[pick] = new_pick
     }
     return(individual)
}
         
         
    
    


    