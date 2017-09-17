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
transformed data {
    // Not put much thought into this prior - just used the Jeffreys prior.
    // One thing to try would be breaking the symmetry somehow by using random
    // perturbation of a weak prior on exposures to try to improve mixing.
    vector[S] kappa = rep_vector(0.5, S);
}
parameters {
    simplex[C] signatures[S];   // matrix of signatures, with simplex constraint
    simplex[S] exposures[G];
    vector<lower=0>[G] multiplier;
}
transformed parameters {
    // Poisson parameters
    // array_to_matrix is defined in common_functions.stan and is not in base Stan
    matrix[G, C] lambda = array_to_matrix(exposures) * array_to_matrix(signatures) .* opps;
    for (g in 1:G) {
        lambda[g] = lambda[g] * multiplier[g];
    }
}
model {
    // The problem with this model is that the probabilities
    // are invariant to perturbations of the columns of
    // exposures and the rows of signatures (the same perturbation)
    // applied to both. This is the old "label switching" problem.
    // Currently nothing is being done about it in the model,
    // and probably avoiding running multiple chains will be fine.
    
    for (s in 1:S) {
        // Priors for signatures
        signatures[s] ~ dirichlet(alpha[s]);
    }

    for (g in 1:G) {
        // Priors for exposures (Jeffreys)
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
