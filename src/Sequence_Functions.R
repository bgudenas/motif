
# Create function to generate random sequences  ---------------------------
library(Biostrings)
library(stringr)
library(flexclust)
library(stringdist)


Rando_DNA_Strings = function(seq_num, width){
#generates a list of random DNA sequence 
    seqs = vector(mode="character", length = seq_num)
    for ( i in 1:seq_num){
    seqs[i] = paste(sample(DNA_ALPHABET[1:4], size=width, replace=TRUE), collapse = "")
    }
    return(seqs)
}


Motif_Implanter = function(seqs, motif, var_num){
    old_motif = motif
    library(stringr)
    seq_length = nchar(seqs[1])
    
     motif_length=nchar(motif)
    for ( i in 1:length(seqs)){
        motif=old_motif
        if ( var_num !=0){
            for (j in 1:var_num){
            var_pos = sample(1:nchar(motif), 1)
            mot_prefix = str_sub(motif, 1, end = var_pos-1)
            new_pos = sample(c("A","T","C","G"), 1)
            mot_suffix = str_sub(motif, var_pos+1, end = nchar(motif))
            motif = str_c(mot_prefix, new_pos, mot_suffix)
            }
            print(stringdist(old_motif, motif))
        }
        
        DNA = seqs[i]
        index = sample(1:(seq_length-motif_length), 1)
        prefix = str_sub(DNA, start = 1, end = index)
        suffix = str_sub(DNA, start = index + motif_length + 1, end = nchar(DNA))
        new_seq = str_c(prefix,motif,suffix)
        seqs[i] = new_seq
    }
return(seqs)
}

Seq2Fasta = function(Seqs){
    sink(paste0("fasta_",nchar(Seqs[1]), ".txt"))
    for ( i in 1:length(Seqs)){
        cat(">")
        cat(paste(i))
        cat("\n")
        cat(Seqs[i])
        cat("\n")
        
    }
    sink()
}
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


PWM_dist = function(individual){
    # pos_set = individual[1:positives]
    # neg_set = individual[(positives+1) : length(individual)]
    
    motif_pos = DNAStringSet(Seqs, start = individual, end = individual + (motif_width-1))
    # motif_neg = DNAStringSet( Seqs[(positives+1) : length(individual) ], start = neg_set, end = neg_set + (motif_width-1))
    PWM_pos = consensusMatrix(motif_pos, as.prob = TRUE)[1:4,]
    IC = 0
    for ( i in 1:nrow(PWM_pos)){
        for (j in 1:ncol(PWM_pos)){
            IC = IC + PWM_pos[i,j]*log2(PWM_pos[i,j] / prior.params[names(prior.params) == rownames(PWM_pos)[i]] )
        }
    }
    return(IC)
}

    
mutation = function(object, parent){
     individual <- as.vector(object@population[parent, ])
     
        pick = sample(1:length(Seqs),1)
     # if (pick <= positives){
         # pos_set = individual[1:positives]
         pos_set = individual[-pick]
         
         motif = DNAStringSet(Seqs[-pick], start = pos_set, end = pos_set + (motif_width-1))
         
         PWM_mat = PWM(motif, prior.params = prior.params)
         
         probs = PWMscoreStartingAt(PWM_mat, Seqs[pick], starting.at = 1:(nchar(Seqs[pick])-motif_width))
         probs[probs < 0 ] = 0
         probs = probs /sum(probs)
         
         new_pick = sample(1:(nchar(Seqs[pick])-motif_width), prob = probs, 1)
         individual[pick] = new_pick
         
     # } else {
     #     neg_set = individual[(positives+1):length(individual)]  ### solution vector
     #     neg_set = neg_set[-(pick-positives)]
     #     
     #     neg_picks = (positives+1):length(Seqs) ## index
     #     neg_picks = neg_picks[-(pick-positives)]
     #     
     #     motif = DNAStringSet(Seqs[neg_picks], start = neg_set, end = neg_set + (motif_width-1))
     #     
     #     PWM_mat = PWM(motif)
     #     
     #     probs = PWMscoreStartingAt(PWM_mat, Seqs[pick], starting.at = 1:(nchar(Seqs[pick])-motif_width))
     #     probs[probs < median(probs)] = 0
     #     probs = probs /sum(probs)
     #     
     #     new_pick = sample(1:(nchar(Seqs[pick])-motif_width), prob = probs, 1)
     #     individual[pick] = new_pick
     # }
     return(individual)
}
         
         
pt_crossover = function(object, parentals){   
    parents <- object@population[parentals, , drop = FALSE]
    children <- matrix(as.double(0), nrow = 2, ncol = ncol(parents))
    fitnessChildren <- rep(NA, 2)
    
    children[1,] = parents[1,]
    children[2,] = parents[2,]
    
    break_point = sample(2: (ncol(children)-2), 1)
    
    children[1, break_point : ncol(parents)] = children[2, break_point : ncol(parents)]
    
    children[2, break_point : ncol(parents)] = parents[1, break_point : ncol(parents)]
    
    out <- list(children = children, fitness = fitnessChildren)
    return(out) 
}
# 
# 
# pwm2ic<-function(pwm) {
#     npos = ncol(pwm)
#     ic = numeric(length=npos)
#     for (i in 1:npos) {
#         ic[i] = 2 + sum(sapply(pwm[, i], function(x) { 
#             print(x)
#             if (x > 0) { x*log2(x) } else { 0 }
#         }))
#     }    
#     return(ic)
# }




pwm2ic <-function(pwm) {
    npos = ncol(pwm)
    ic = numeric(length=npos)
    for (j in 1:npos) {
        for (i in 1:nrow(pwm)) {
            
            if (pwm[i,j] > 0){ 
                prior.prob = prior.params[names(prior.params) == names(pwm[i,j])]
                
                ic[j] = ic[j] + pwm[i,j] * log2(pwm[i, j]/ prior.prob)  
            }
        }
    }
            
    return(sum(ic))
}



dinuc_shuffle = function(sequence){
    shuff=""
    while (nchar(sequence) > 1){
        index = sample(1:(nchar(sequence)-1), 1)
        di_nuc = str_sub(sequence, start = index, end = index+1)
        shuff = str_c(shuff, di_nuc)
        
        prefix = str_sub(string = sequence, start = 1, end = (index-1))
        suffix = str_sub(string = sequence, start = (index+2), end = nchar(sequence))
        sequence = str_c(prefix, suffix)
        # print(nchar(sequence))
    }
    return(shuff)
}


pwm_2_seq = function(pwm){
    seq=""
    
    for ( i in 1:ncol(pwm)){
        let = names(which.max(pwm[ ,i]))
        if (max(pwm[ ,i]) < 0.5){ let = "N"}
        seq = str_c(seq,let)
    }
    return(seq)
}
# pwm_2_seq(pwm)       
#     
# string_dist()
#     


    