data {
    int<lower=1> C;  // number of mutation categories
    int<lower=1> S;  // number of mutational signatures
    int<lower=1> G;  // number of genomes
    matrix[S, C] signatures;  // matrix of signatures (columns) to be fitted
    int counts[G, C];         // data = counts per category (columns) per genome sample (rows)
    vector<lower=0>[S] alpha; // prior on exposures (i.e. mixing proportions of signatures)
    int opps[G, C];       // Matrix of opportunities
}
parameters {
    simplex[S] exposures[G];
    real<lower=0> multiplier[G];
}
model {
    {
        // Set parameters for likelihood
        matrix[G, S] scaled_exposures;
        matrix[G, C] lambda; // poisson parameters
        for (i in 1:G) {
            scaled_exposures[i] = to_row_vector(exposures[i]) * multiplier[i];
        }
        lambda = (scaled_exposures * signatures) .* to_matrix(opps);
        
        // Priors
        for (i in 1:G) {
            exposures[i] ~ dirichlet(alpha);
            multiplier ~ cauchy(0, 1);
        }
        
        // Likelihood
        for (i in 1:G) {
            counts[i] ~ poisson(lambda[i]);
        }
    }
}
generated quantities {
    vector[G*C] log_lik;
    {
        matrix[G, S] scaled_exposures;
        matrix[G, C] lambda; // poisson parameters

        for (i in 1:G) {
            scaled_exposures[i] = to_row_vector(exposures[i]) * multiplier[i];
        }

        lambda = (scaled_exposures * signatures) .* to_matrix(opps);
        for (i in 1:G) {
            for (j in 1:C) {
                log_lik[(i-1)*C + j] = poisson_lpmf(counts[i, j] | lambda[i, j]);
            }
        }
    }
}
