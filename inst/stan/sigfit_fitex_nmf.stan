functions {
#include /include/common_functions.stan
}
data {
    int<lower=1> C;            // number of categories
    int<lower=1> S;            // number of fixed signatures
    int<lower=1> G;            // number of genomes
    int<lower=1> N;            // number of extra signatures
    matrix[S, C] fixed_sigs;   // matrix of signatures (rows) by categories (columns)
    int<lower=0> counts[G, C]; // matrix of counts per genome (rows) in each category (columns)
    matrix[N, C] alpha;        // priors for extra signatures
    real<lower=0> kappa;       // prior for exposures
}
transformed data {
    int T = S + N;   // total number of signatures, including extra signatures
    vector[T] kappa_vec = rep_vector(kappa, T);
}
parameters {
    simplex[C] extra_sigs[N];   // additional signatures to extract
    simplex[T] exposures[G]; // includes exposures for extra_sigs
}
transformed parameters {
    matrix[T, C] signatures;
    matrix<lower=0, upper=1>[G, C] probs;

    signatures = append_row(fixed_sigs, array_to_matrix(extra_sigs));
    probs = array_to_matrix(exposures) * signatures;
}
model {
    // Priors for extra signatures
    for (n in 1:N) {
        extra_sigs[n] ~ dirichlet(alpha[n]');
    }

    for (g in 1:G) {
        // Priors for exposures_raw are uniform dirichlet
        exposures[g] ~ dirichlet(kappa_vec);

        // Likelihood
        counts[g] ~ multinomial(to_vector(probs[g]));
    }
}
