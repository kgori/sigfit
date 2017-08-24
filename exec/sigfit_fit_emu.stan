functions {
    #include "common_functions.stan"
}
data {
    int<lower=1> C;  // number of mutation categories
    int<lower=1> S;  // number of mutational signatures
    int<lower=1> G;  // number of genomes
    matrix[S, C] signatures;  // matrix of signatures (columns) to be fitted
    int<lower=0> counts[G, C];         // data = counts per category (columns) per genome sample (rows)
    vector<lower=0>[S] alpha; // prior on exposures (i.e. mixing proportions of signatures)
    matrix[G, C] opps;        // Matrix of opportunities
}
parameters {
    simplex[S] exposures[G];
    real<lower=0> multiplier[G];
}
transformed parameters {
    // poisson parameters
    matrix[G, C] lambda = array_to_matrix(exposures) * signatures .* opps;
    for (g in 1:G) {
        lambda[g] = lambda[g] * multiplier[g];
    }
}
model {
    for (i in 1:G) {
        // Priors
        exposures[i] ~ dirichlet(alpha);
        multiplier ~ cauchy(0, 1);
        
        // Likelihood
        counts[i] ~ poisson(lambda[i]);
    }
}
generated quantities {
    vector[G*C] log_lik;
    real bic;
    for (i in 1:G) {
        for (j in 1:C) {
            log_lik[(i-1)*C + j] = poisson_lpmf(counts[i, j] | lambda[i, j]);
        }
    }
    bic = 2 * sum(log_lik) - log(G) * (G*S);
}
