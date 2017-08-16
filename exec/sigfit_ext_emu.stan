data {
    int<lower=1> C;  // number of mutation categories [uses index c]
    int<lower=1> S;  // number of mutational signatures [uses index s]
    int<lower=1> G;  // number of genomes [uses index g]
    int counts[G, C];  // Matrix of mutation counts per sample (rows) per category
    matrix[G, C] opps; // Matrix of opportunities per sample (rows) per category
}
transformed data {
    vector[C] alpha;
    for (c in 1:C) {
        alpha[c] = 1;
    }
}
parameters {
    simplex[C] signatures[S];           // Matrix of signatures, with simplex constraint
    matrix<lower=0>[G, S] exposures;
}
transformed parameters {
    // Precomputation
    matrix[G, C] lambda;    // elementwise product of sig_expo and opps
    {
        matrix[S, C] signatures_mat;    // signatures recast as matrix
    
        for (s in 1:S) {
            for (c in 1:C) {
                signatures_mat[s, c] = signatures[s, c];
            }
        }
        lambda = exposures * signatures_mat .* opps;
    }
}
model {
    // Priors for exposures
    for (g in 1:G) {
        exposures[g] ~ cauchy(0, 2.5); // need something positive continuous
    }
    
    // Priors for signatures
    for (s in 1:S) {
        signatures[s] ~ dirichlet(alpha);
    }
    
    // Likelihood
    for (g in 1:G) {
        counts[g] ~ poisson(lambda[g]);
    }
}
generated quantities {
    vector[G*C] log_lik;
    real bic;
    
    // Compute log_lik
    for (i in 1:G) {
        for (j in 1:C) {
            // TODO: vectorise this?
            log_lik[(i-1)*C + j] = poisson_lpmf(counts[i, j] | lambda[i, j]);
        }
    }
    
    // Compute bic with (G*S + S*(C-1)) free parameters
    bic = 2 * sum(log_lik) - log(G) * (G*S + S*(C-1));
}
