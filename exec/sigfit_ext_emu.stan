functions {
    #include "common_functions.stan"
}
data {
    int<lower=1> C;             // number of mutation categories [uses index c]
    int<lower=1> S;             // number of mutational signatures [uses index s]
    int<lower=1> G;             // number of genomes [uses index g]
    int<lower=0> counts[G, C];  // matrix of mutation counts per sample (rows) per category
    matrix[G, C] opps;          // matrix of opportunities per sample (rows) per category
    matrix[S, C] alpha;         // prior for signatures
}
parameters {
    simplex[C] signatures[S];   // matrix of signatures, with simplex constraint
    matrix<lower=0>[G, S] activities;
}
transformed parameters {
    // Poisson parameters
    // array_to_matrix is defined in common_functions.stan and is not in base Stan
    matrix[G, C] expected_counts = activities * array_to_matrix(signatures) .* opps;
}
model {
    // Label Switching -  the probabilities
    // are invariant to perturbations of the columns of
    // exposures and the rows of signatures (the same perturbation
    // applied to both). Avoiding running multiple chains.
    
    for (s in 1:S) {
        // Priors for signatures
        signatures[s] ~ dirichlet(alpha[s]');
    }

    for (g in 1:G) {
        // Priors for activities
        activities[g] ~ cauchy(0, 1);

        // Likelihood
        counts[g] ~ poisson(expected_counts[g]);
    }
}
generated quantities {
    matrix[G, S] exposures;
    for (g in 1:G) {
        exposures[g] = scale_row_to_sum_1(activities[g]);
    }
}
