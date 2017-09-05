functions {
    #include "common_functions.stan"
}
data {
    int<lower=1> C;  // number of mutation categories [uses index c]
    int<lower=1> S;  // number of mutational signatures [uses index s]
    int<lower=1> G;  // number of genomes [uses index g]
    int<lower=0> counts[G, C];  // Matrix of mutation counts per sample (rows) per category
    matrix[G, C] opps; // Matrix of opportunities per sample (rows) per category
}
transformed data {
    // uniform priors for Dirichlets
    vector[C] alpha = rep_vector(1, C);
    vector[S] kappa = rep_vector(1, S);
}
parameters {
    simplex[C] signatures[S]; // Matrix of signatures, with simplex constraint
    simplex[S] exposures[G];
    vector<lower=0>[G] multiplier;
}
transformed parameters {
    // Poisson parameters
    matrix[G, C] lambda = array_to_matrix(exposures) * array_to_matrix(signatures) .* opps;
    for (g in 1:G) {
        lambda[g] = lambda[g] * multiplier[g];
    }
}
model {
    // Priors for signatures
    for (s in 1:S) {
        signatures[s] ~ dirichlet(alpha);
    }

    for (g in 1:G) {
        // Priors for exposures
        exposures[g] ~ dirichlet(kappa);
        multiplier[g] ~ cauchy(0, 1);
        
        // Likelihood
        counts[g] ~ poisson(lambda[g]);
    }
}
generated quantities {
    vector[G] log_lik;
    real bic;
    
    // Compute log_lik
    for (g in 1:G) {
        log_lik[g] = poisson_lpmf(counts[g] | lambda[g]);
    }
    
    // Compute bic with (G*S + S*(C-1)) free parameters
    bic = 2 * sum(log_lik) - log(G) * (G*S + S*(C-1));
}
