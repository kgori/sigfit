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
    matrix<lower=0>[G,S] exposures_raw;
}
transformed parameters {
    // Poisson parameters
    // array_to_matrix is defined in common_functions.stan and is not in base Stan
    matrix[G, C] lambda = exposures_raw * array_to_matrix(signatures) .* opps;
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
        // Priors for exposures
        exposures_raw[g] ~ cauchy(0, 1);

        // Likelihood
        counts[g] ~ poisson(lambda[g]);
    }
}
generated quantities {
    matrix[G, S] exposures;
    vector[G] log_lik;
    
    // reconstructed sample counts - do counts estimated by the model look like the input data?
    matrix[S, C] reconstruction[G];
    
    for (g in 1:G) {
        exposures[g] = scale_row_to_sum_1(exposures_raw[g]);
        
        log_lik[g] = poisson_lpmf(counts[g] | lambda[g]);
        
        for (s in 1:S) {
            reconstruction[g][s] = (exposures_raw[g, s] * signatures[s])' .* opps[g];
        }
    }
}
