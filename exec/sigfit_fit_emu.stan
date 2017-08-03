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
    simplex[S] exposure_probs[G];
    real<lower=0> multiplier;
}
transformed parameters {
    matrix[G, S] exposures;
    for (i in 1:G) {
        exposures[i] = to_row_vector(exposure_probs[i]) * multiplier;
    }
}
model {
    {
        matrix[G, C] lambda; // poisson parameters
        lambda = (exposures * signatures) .* to_matrix(opps);
        
        for (i in 1:G) {
            exposure_probs[i] ~ dirichlet(alpha);
        }
        multiplier ~ cauchy(0, 1);
        
        // Likelihood
        for (i in 1:G) {
            counts[i] ~ poisson(lambda[i]);
        }
    }
}
