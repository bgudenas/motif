
library(gtools)

seeder = function(object) {
    pop_size = object@popSize
    
    max_sizes = unlist(lapply(Seqs, function(x) nchar(x)))
    seeds = matrix( ncol = 4*motif_width, nrow = pop_size)
    for ( i in 1:pop_size){
        seeds[i,] = as.vector(t(rdirichlet(motif_width, c(1,1,1,1))))
    }
    return(seeds)
}



PWM_fitness = function(individual){
    pwm = matrix(individual, nrow=4, ncol = motif_width)
    rownames(pwm) = c("A", "C", "G", "T")
    
    IC = sum(pwm2ic(pwm))
    pwm = pwm+ 0.0001
    
    negatives = (length(Seqs)-positives)
    rand_seqs = sample(negatives :length(Seqs), size = 0.2*negatives)
    scores = vector(mode = "numeric", length = length(rand_seqs))
    for ( rand in 1:length(rand_seqs)) {
        i= rand_seqs[rand]
        scores[rand] = countPWM(pwm, Seqs[i], min.score = "70%")
    }
    
    # pwm = pwm + 0.000001
    # pwm = log(pwm)
    # 
    # scores = vector(mode = "numeric", length = length(Seqs))
    # 
    # best = sum(apply(pwm, 2, max))
    # if (best < 0){
    #     threshold = best + (best*0.33)
    # } else threshold = best*0.33
    # 
    # for ( i in 1:length(Seqs)){
    #     scores_seq  = PWMscoreStartingAt(pwm, Seqs[[i]], starting.at = 1 : (nchar(Seqs[[i]]) - motif_width ) )
    #     scores_seq[scores_seq < threshold] = 0
    #     scores[i] = sum(scores_seq)
    # }
    # 
    # fitness = abs(  sum(scores[1:positives]) - sum(scores[(positives+1):length(scores)]) )
    
    return(IC/(sum(scores)+1))
}





PWM_mutation = function(object, parent){
    individual <- as.vector(object@population[parent, ])
    pwm = matrix(individual, nrow=4, ncol = motif_width)
    rownames(pwm) = c("A", "C", "G", "T")
    
    # maxes = apply(pwm, 2, max)
    # 
    # maxes = 1-maxes
    # probs = maxes/sum(maxes)
    # 
    # mut_chose = sample(1:ncol(pwm), 1, prob = probs)
    # 
    # pwm[ , mut_chose] = pwm[, mut_chose][sample(1:4)]
    
    pwm = pwm + 0.0001
    pwm = log(pwm)

    # scores = vector(mode = "numeric", length = length(Seqs))
# 
#     best = sum(apply(pwm, 2, max))
#     if (best < 0){
#         threshold = best + (best*0.33)
#     } else threshold = best*0.33
    rand_seqs = sample(1:positives, size = 0.3*positives)
    scores = vector(mode = "numeric", length = length(rand_seqs))
    for ( rand in 1:length(rand_seqs)) {
        i= rand_seqs[rand]
        scores_seq  = PWMscoreStartingAt(pwm, Seqs[[i]], starting.at = 1 : (nchar(Seqs[[i]]) - motif_width ) )
        # if (i <= positives){
        #     scores[i] = which.max(scores_seq)
        # }else { scores[i] = which.max(scores_seq)
        # }
        
        scores[rand] = which.max(scores_seq)
    }
    mut_pos = DNAStringSet(Seqs[rand_seqs], start = scores, end = scores + (motif_width-1))
    pwm = consensusMatrix(mut_pos, as.prob = TRUE)[1:4, ]
    
    return(as.vector(pwm))
}

PWM_crossover = function(object, parentals){   
    parents <- object@population[parentals, , drop = FALSE]
    fitnessChildren <- rep(NA, 2)
    
    mut_num = sample(2:(motif_width-2), 1)  ## chose num of cols to mutate (2:motif_width-2)
    
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
    
    
    
                       
                       