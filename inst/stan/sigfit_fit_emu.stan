functions {
#include /include/common_functions.stan
}
data {
    int<lower=1> C;             // number of mutation categories
    int<lower=1> S;             // number of mutational signatures
    int<lower=1> G;             // number of genomes
    int<lower=0> counts[G, C];  // observed mutation counts (genome per row)
    matrix[S, C] signatures;    // signatures to fit (signature per row)
    vector<lower=0>[S] kappa;   // prior on exposures (mixing proportions)
    matrix[G, C] opps;          // mutational opportunities (genome per row)
}
parameters {
    simplex[S] exposures[G];      // signature exposures (genome per row)
    real<lower=0> multiplier[G];  // exposure multipliers
}
transformed parameters {
    matrix<lower=0>[G, S] activities;  // signature activities (genome per row)
    matrix[G, C] expected_counts;      // Poisson parameters (genome per row)
    for (g in 1:G) {
        activities[g] = exposures[g]' * multiplier[g];
    }
    expected_counts = activities * signatures .* opps;
}
model {
    multiplier ~ cauchy(0, 1);                    // multiplier prior
    for (g in 1:G) {
        exposures[g] ~ dirichlet(kappa);          // exposure priors
        counts[g] ~ poisson(expected_counts[g]);  // Poisson likelihood
    }
}
