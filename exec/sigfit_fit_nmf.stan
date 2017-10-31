functions {
    #include "common_functions.stan"
}
data {
    int<lower=1> C;             // number of mutation categories
    int<lower=1> S;             // number of mutational signatures
    int<lower=1> G;             // number of genomes
    matrix[S, C] signatures;    // matrix of signatures (rows) to be fitted
    int<lower=0> counts[G, C];  // matrix of counts per category (columns) per genome sample (rows)
    vector<lower=0>[S] kappa;   // prior on exposures (i.e. mixing proportions of signatures)
}
parameters {
    simplex[S] exposures[G];
}
transformed parameters {
    // array_to_matrix is defined in common_functions.stan and is not in base Stan
    matrix<lower=0, upper=1>[G, C] probs = array_to_matrix(exposures) * signatures;
}
model {
    for (g in 1:G) {
        // Priors
        exposures[g] ~ dirichlet(kappa);
        
        // Likelihood
        counts[g] ~ multinomial(probs[g]');
    }
}
