data {
    int<lower=1> C;  // number of mutation categories
    int<lower=1> S;  // number of mutational signatures
    int<lower=1> G;  // number of genomes
    matrix[S, C] signatures;  // matrix of signatures (columns) to be fitted
    int counts[G, C];         // data = counts per category (columns) per genome sample (rows)
    vector<lower=0>[S] alpha; // prior on exposures (i.e. mixing proportions of signatures)
    int opps[G, C];           // Matrix of opportunities
}
parameters {
    simplex[S] exposure_probs[G];
    real<lower=0> multiplier;
}
transformed parameters {
    matrix[G, S] exposures;
    for (g in 1:G) {
        exposures[g] = to_row_vector(exposure_probs[g]) * multiplier;
    }
}
model {
    {
        matrix[G, C] lambda; // poisson parameters
        lambda = (exposures * signatures) .* to_matrix(opps);
        
        for (g in 1:G) {
            exposure_probs[g] ~ dirichlet(alpha);
        }
        multiplier ~ cauchy(0, 1);
        
        // Likelihood
        for (g in 1:G) {
            counts[g] ~ poisson(lambda[g]);
        }
    }
}
generated quantities {
    vector[G] log_lik;
    for (g in 1:G) {
        log_lik[g] = poisson_lpmf(counts[g] | lambda[g]);
    }
}
