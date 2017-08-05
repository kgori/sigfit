functions {
    #include "common_functions.stan"
}
data {
    int<lower=1> C;  // number of mutation categories
    int<lower=1> S;  // number of mutational signatures
    int<lower=1> G;  // number of genomes
    matrix[S, C] signatures;  // matrix of signatures (rows) to be fitted
    int counts[G, C];         // data = counts per category (columns) per genome sample (rows)
    vector<lower=0>[S] alpha; // prior on exposures (i.e. mixing proportions of signatures)
}
parameters {
    simplex[S] exposures[G];
}
transformed parameters {
    matrix[G, S] exposures_mat;
    matrix<lower=0, upper=1>[G, C] probs;
    for (i in 1:G) {
        for (j in 1:S) {
            exposures_mat[i, j] = exposures[i, j];
        }
    }
    probs = exposures_mat * signatures;
}
model {
    for (i in 1:G) {
        exposures[i] ~ dirichlet(alpha);
        counts[i] ~ multinomial(probs[i]);
    }
}
generated quantities {
    vector[G] log_lik;
    for (i in 1:G) {
        log_lik[i] = multinomial_lpmf(counts[i] | probs[i]);
    }
}
