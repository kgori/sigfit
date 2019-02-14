functions {
    #include "common_functions.stan"
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
    matrix<lower=0>[G, C] phi;    // negative binomial overdispersion
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
    multiplier ~ cauchy(0, 1);            // multiplier priors
    for (g in 1:G) {
        exposures[g] ~ dirichlet(kappa);  // exposure priors
        phi[g] ~ cauchy(0, 2.5);          // overdispersion priors
        for (c in 1:C) {
            counts[g, c] ~ neg_binomial_2(expected_counts[g, c], phi[g, c]);  // neg binomial likelihood
        }
    }
}
