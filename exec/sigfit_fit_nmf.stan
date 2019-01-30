functions {
    #include "common_functions.stan"
}
data {
    int<lower=1> C;             // number of mutation categories
    int<lower=1> S;             // number of mutational signatures
    int<lower=1> G;             // number of genomes
    int<lower=0> counts[G, C];  // observed mutation counts (row per genome)
    matrix[S, C] signatures;    // signatures to fit (row per signature)
    vector<lower=0>[S] kappa;   // prior on exposures (mixing proportions)
}
parameters {
    simplex[S] exposures[G];  // signature exposures (row per genome)
}
transformed parameters {
    // array_to_matrix is defined in common_functions.stan and is not in base Stan
    // Mutation type probabilities (row per genome)
    matrix<lower=0, upper=1>[G, C] probs = array_to_matrix(exposures) * signatures;
}
model {
    for (g in 1:G) {
        exposures[g] ~ dirichlet(kappa);     // exposure priors
        counts[g] ~ multinomial(probs[g]');  // multinomial likelihood
    }
}
