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
    matrix[G, C] opps;         // matrix of opportunities per sample (rows) per category
    matrix[N, C] alpha;        // priors for extra signatures
}
transformed data {
    int T = S + N;   // total number of signatures, including extra signatures
}
parameters {
    simplex[C] extra_sigs[N];             // additional signatures to extract
    matrix<lower=0>[G, T] activities;  // includes exposures for extra_sigs
}
transformed parameters {
    // Full signatures matrix
    matrix[T, C] signatures = append_row(fixed_sigs, array_to_matrix(extra_sigs));
    
    // Poisson parameters
    matrix[G, C] expected_counts = activities * signatures .* opps;
}
model {
    // Priors for extra signatures
    for (n in 1:N) {
        extra_sigs[n] ~ dirichlet(alpha[n]');
    }

    for (g in 1:G) {
        // Priors for activities
        activities[g] ~ cauchy(0, 1);
        
        // Likelihood
        counts[g] ~ poisson(expected_counts[g]);
    }
}
generated quantities {
    matrix[G, T] exposures;
    for (g in 1:G) {
        exposures[g] = scale_row_to_sum_1(activities[g]);
    }
}
