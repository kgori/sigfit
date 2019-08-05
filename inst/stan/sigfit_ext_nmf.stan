functions{
#include /include/common_functions.stan
}
data {
    int<lower=1> C;             // number of mutation categories
    int<lower=1> S;             // number of mutational signatures
    int<lower=1> G;             // number of genomes
    int<lower=0> counts[G, C];  // matrix of counts per category (columns) per genome sample (rows)
    matrix[S, C] alpha;         // prior for signatures
    real<lower=0> kappa;        // prior for exposures
}
transformed data {
    vector<lower=0>[S] kappa_vec = rep_vector(kappa, S);
}
parameters {
    simplex[S] exposures[G];
    simplex[C] signatures[S];
}
transformed parameters {
    matrix<lower=0, upper=1>[G, C] probs;
    // array_to_matrix is defined in common_functions.stan and is not in base Stan
    probs = array_to_matrix(exposures) * array_to_matrix(signatures);
}
model {
    for (s in 1:S) {
        // Priors for signatures
        signatures[s] ~ dirichlet(alpha[s]');
    }

    for (g in 1:G) {
        // Priors for exposures uniform dirichlet
        exposures[g] ~ dirichlet(kappa_vec);

        // Likelihood
        counts[g] ~ multinomial(probs[g]');
    }
}
