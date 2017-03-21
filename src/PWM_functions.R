library(gtools)
library(Biostrings)


priors = letterFrequency(DNAStringSet(Seqs), letters = DNA_BASES)
prior.params = colSums(priors)/sum(priors)
 
seeder = function(object) {
    pop_size = object@popSize
    
    # max_sizes = unlist(lapply(Seqs, function(x) nchar(x)))
    seeds = matrix( ncol = 4*motif_width, nrow = pop_size)
    for ( i in 1:pop_size){
        seeds[i,] = as.vector(t(rdirichlet(motif_width, c(1,1,1,1))))
    }
    return(seeds)
}


PWM_fitness = function(individual){
    pwm = matrix(individual, nrow=4, ncol = motif_width)
    rownames(pwm) = c("A", "C", "G", "T")
    
    IC = pwm2ic(pwm)
    pwm = pwm + 0.0001
    
    negatives = (length(Seqs)-positives)+1
    rand_seqs = sample(negatives:length(Seqs), size = sqrt(length(negatives))) #0.05*negatives
    scores = vector(mode = "numeric", length = length(rand_seqs))
    # neg_scores = vector(mode = "numeric", length = length(rand_seqs))
    for ( rand in 1:length(rand_seqs)) {
        i = rand_seqs[rand]
        norm_score = PWMscoreStartingAt(pwm, Seqs[i], 1:(nchar(Seqs[i])-motif_width))
        # rev = reverseComplement(DNAString(Seqs[i]))
        # rev_score = countPWM(pwm, rev, min.score = "60%")
       scores[rand] = mean(norm_score)/maxScore(pwm)
       
    }

    return(IC/(sum(scores)))
}




rando_mutation = function(object, parent){
    individual <- as.vector(object@population[parent, ])
    pwm = matrix(individual, nrow=4, ncol = motif_width)
    rownames(pwm) = c("A", "C", "G", "T")
    
    mut_pos = sample(1:ncol(pwm), 1)
    pwm[ ,mut_pos] = as.vector(t(rdirichlet(1, c(1,1,1,1))))
    return(as.vector(pwm))
}
    


PWM_mutation = function(object, parent){
    individual <- as.vector(object@population[parent, ])
    pwm = matrix(individual, nrow=4, ncol = motif_width)
    rownames(pwm) = c("A", "C", "G", "T")
    
    pwm = pwm + 0.0001
    pwm = log2(pwm)

    rand_seqs = sample(1:positives, size = sqrt(positives)) #0.01*positives
    rSeqs = Seqs[rand_seqs]
    scores = vector(mode = "numeric", length = length(rand_seqs))
    for ( rand in 1:length(rand_seqs)) {
        # i= rand_seqs[rand]
        scores_seq  = PWMscoreStartingAt(pwm, rSeqs[rand], starting.at = 1 : (nchar(rSeqs[rand]) - motif_width ) )
        
        scores_rev = PWMscoreStartingAt(pwm, reverseComplement(DNAString(rSeqs[rand]), starting.at = 1 : (nchar(rSeqs[rand])) - motif_width ) )

        if ( max(scores_rev) > max(scores_seq) ){
            max_score = which.max(scores_rev)
            rSeqs[rand] = as.character(reverseComplement(DNAString(rSeqs[rand])))

        } else {max_score = which.max(scores_seq)}

        scores[rand] = max_score
    }
    mut_pos = DNAStringSet(rSeqs, start = scores, end = scores + (motif_width-1))
    pwm = consensusMatrix(mut_pos, as.prob = TRUE)[1:4, ]
    
    return(as.vector(pwm))
}

PWM_crossover = function(object, parentals){   
    parents <- object@population[parentals, , drop = FALSE]
    fitnessChildren <- rep(NA, 2)
    
    mut_num = sample(1:(motif_width-1), 1)  ## chose num of cols to mutate (2:motif_width-2)
    
    mut_pos = sample(1:motif_width, mut_num) # randomly select mut_num of cols to swap
    
    child1 = matrix(parents[1,], nrow=4, ncol = motif_width)
    child2 =  matrix(parents[2,], nrow=4, ncol = motif_width)
    
    child1_tmp = child1[ , mut_pos]
    child1[ , mut_pos] = child2[ , mut_pos]
    
    child2[ , mut_pos] = child1_tmp
    
    
    children = rbind(as.vector(child1), as.vector(child2))
    
    out = list(children = children, fitness = fitnessChildren)
    return(out) 
}
    
    
    
                       
                       