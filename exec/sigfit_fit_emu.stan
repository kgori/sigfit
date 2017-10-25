functions {
    #include "common_functions.stan"
}
data {
    int<lower=1> C;            // number of mutation categories
    int<lower=1> S;            // number of mutational signatures
    int<lower=1> G;            // number of genomes
    matrix[S, C] signatures;   // matrix of signatures (columns) to be fitted
    int<lower=0> counts[G, C]; // data = counts per category (columns) per genome sample (rows)
    vector<lower=0>[S] alpha;  // prior on exposures (i.e. mixing proportions of signatures)
    matrix[G, C] opps;         // matrix of opportunities
}
parameters {
    simplex[S] exposures[G];
    real<lower=0> multiplier[G];
}
transformed parameters {
    matrix<lower=0>[G, S] exposures_raw;
    matrix[G, C] lambda;  // Poisson parameters
    
    for (g in 1:G) {
        exposures_raw[g] = exposures[g]' * multiplier[g];
    }

    // Poisson parameters
    lambda = exposures_raw * signatures .* opps;
}
model {
    for (i in 1:G) {
        // Priors
        exposures[i] ~ dirichlet(alpha);
        multiplier ~ cauchy(0, 1);
        
        // Likelihood
        counts[i] ~ poisson(lambda[i]);
    }
}
