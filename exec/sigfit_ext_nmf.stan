functions{
    #include "common_functions.stan"
}
data {
    int<lower=1> C;             // number of mutation categories
    int<lower=1> S;             // number of mutational signatures
    int<lower=1> G;             // number of genomes
    int<lower=0> counts[G, C];  // matrix of counts per category (columns) per genome sample (rows)
    matrix[S, C] alpha;         // prior for signatures
}
parameters {
    vector[S] exposures_raw[G];
    simplex[C] signatures[S];
}
transformed parameters {
    matrix<lower=0, upper=1>[G, S] exposures;
    matrix<lower=0>[G, C] probs;
    
    // Convert offset parameter(s) to +ve, sum-to-one vectors
    for (g in 1:G) {
        exposures[g] = softmax(exposures_raw[g])';
    }
    
    // array_to_matrix is defined in common_functions.stan and is not in base Stan
    probs = exposures * array_to_matrix(signatures); 
}
model {
    for (s in 1:S) {
        // Priors for signatures
        signatures[s] ~ dirichlet(alpha[s]');
    }
    
    for (g in 1:G) {
        // Priors for exposures_raw are standard normal
        exposures_raw[g] ~ normal(0, 1);
        
        // Likelihood
        counts[g] ~ multinomial(probs[g]');
    }
}
