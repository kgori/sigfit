functions {
    #include "common_functions.stan"
}
data {
    int<lower=1> C;            // number of categories
    int<lower=1> S;            // number of fixed signatures
    int<lower=1> G;            // number of genomes
    int<lower=1> N;            // number of extra signatures
    matrix[S, C] fixed_sigs;   // matrix of signatures (rows) by categories (columns)
    int<lower=0> counts[G, C]; // matrix of counts per genome (rows) in each category (columns)
    matrix[N, C] alpha;        // priors for extra signatures
}
transformed data {
    int T = S + N;   // total number of signatures, including extra signatures
}
parameters {
    simplex[C] extra_sigs[N];  // additional signatures to extract
    vector[T] exposures_raw[G]; // includes exposures for extra_sigs
}
transformed parameters {
    matrix[G, T] exposures;
    matrix[T, C] signatures;
    matrix<lower=0, upper=1>[G, C] probs;
    
    for (g in 1:G) {
        exposures[g] = softmax(exposures_raw[g])';
    }
    
    signatures = append_row(fixed_sigs, array_to_matrix(extra_sigs));
    probs = exposures * signatures;
}
model {
    // Priors for extra signatures
    for (n in 1:N) {
        extra_sigs[n] ~ dirichlet(alpha[n]');
    }
    
    for (g in 1:G) {
        // Priors for exposures_raw
        exposures_raw[g] ~ normal(0, 1);
        
        // Likelihood
        counts[g] ~ multinomial(to_vector(probs[g]));
    }
}
