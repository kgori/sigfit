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
    matrix[G, C] sig_expo_opps;    // elementwise product of sig_expo and opps
    {
        matrix[S, C] signatures_mat;    // signatures recast as matrix
    
        for (s in 1:S) {
            for (c in 1:C) {
                signatures_mat[s, c] = signatures[s, c];
            }
        }
        sig_expo_opps = exposures * signatures_mat .* opps;
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
        counts[g] ~ poisson(sig_expo_opps[g]);
    }
}
generated quantities {
    vector[G] log_lik;
    for (g in 1:G) {
        log_lik[g] = poisson_lpmf(counts[g] | sig_expo_opps[g]);
    }
    //matrix[G, C] log_lik;
    //for (g in 1:G) {
    //    for (c in 1:C) {
    //        log_lik[g, c] = poisson_lpmf(counts[g, c] | sig_expo_opps[g, c]);
    //    }
    //}
}
