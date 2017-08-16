data {
    int<lower=1> C;  // number of mutation categories [uses index c]
    int<lower=1> S;  // number of mutational signatures [uses index s]
    int<lower=1> G;  // number of genomes [uses index g]
    int counts[G, C];  // Matrix of mutation counts per sample (rows) per category
    matrix[G, C] opps; // Matrix of opportunities per sample (rows) per category
}
transformed data {
    // uniform priors for Dirichlets
    vector[C] alpha;
    vector[S] kappa;
    
    for (c in 1:C) {
        alpha[c] = 1;
    }
    
    for (s in 1:S) {
        kappa[s] = 1;
    }
}
parameters {
    simplex[C] signatures[S]; // Matrix of signatures, with simplex constraint
    simplex[S] exposures[G];
    vector<lower=0>[G] multiplier;
}
transformed parameters {
    // Precomputation
    matrix[G, C] lambda; // Poisson parameters
    {                 
        matrix[S, C] signatures_mat;    // signatures recast as matrix
        matrix[G, S] exposures_mat;    // exposures recast as matrix
        
        for (s in 1:S) {
            for (c in 1:C) {
                signatures_mat[s, c] = signatures[s, c];
            }
        }
        
        for (g in 1:G) {
            for (s in 1:S) {
                exposures_mat[g, s] = exposures[g, s] * multiplier[g];
            }
        }
        lambda = exposures_mat * signatures_mat .* opps;
    }
}
model {
    // Priors for exposures
    for (g in 1:G) {
        exposures[g] ~ dirichlet(kappa);
        multiplier[g] ~ cauchy(0, 1);
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
