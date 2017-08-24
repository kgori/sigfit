functions{
    #include "common_functions.stan"
}
data {
    int<lower=1> C;  // number of mutation categories
    int<lower=1> S;  // number of mutational signatures
    int<lower=1> G;  // number of genomes
    int<lower=0> counts[G, C];  // data = counts per category (columns) per genome sample (rows)
    real<lower=0> exposures_prior_val;
}
transformed data {
    // Not put much thought into these priors - just used the Jeffreys prior.
    // One thing to try would be breaking the symmetry somehow by using random
    // perturbation of a weak prior on exposures to try to improve mixing.
    // Another idea is to allow use of (e.g.) cosmic signatures as priors
    // on the signatures.
    vector[S] exposures_prior = rep_vector(exposures_prior_val, S);
    vector[C] signatures_prior = rep_vector(0.5, C);
}
parameters {
    simplex[S] exposures[G];
    simplex[C] signatures[S];
}
transformed parameters {
    // array_to_matrix is defined in common_functions.stan and is not in base Stan
    matrix<lower=0>[G, C] probs = array_to_matrix(exposures) * array_to_matrix(signatures); 
}
model {
    // The problem with this model is that the probabilities
    // are invariant to perturbations of the columns of
    // exposures and the rows of signatures (the same perturbation)
    // applied to both. This is the old "label switching" problem.
    // Currently I'm not doing anything about it in the model,
    // and probably avoiding running multiple chains will be fine.
    for (s in 1:S) {
        signatures[s] ~ dirichlet(signatures_prior);
    }
    
    for (g in 1:G) {
        exposures[g] ~ dirichlet(exposures_prior);
        counts[g] ~ multinomial(probs[g]');
    }
}
generated quantities {
    vector[G] log_lik;
    real bic;
    
    // Compute log_lik
    for (g in 1:G) {
        log_lik[g] = multinomial_lpmf(counts[g] | probs[g]');
    }
    
    // Compute bic
    bic = 2 * sum(log_lik) - log(G) * (G*(S-1) + S*(C-1));
}
