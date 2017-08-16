data {
    int<lower=1> C;  // number of mutation categories
    int<lower=1> S;  // number of mutational signatures
    int<lower=1> G;  // number of genomes
    int<lower=0> counts[G, C];  // data = counts per category (columns) per genome sample (rows)
    real<lower=0> exposures_prior_val;
    // real<lower=0> signatures_prior[S, C];
}
transformed data {
    vector[S] exposures_prior;
    vector[C] signatures_prior;

    // Not put much thought into these priors - just used the Jeffreys prior.
    // One thing to try would be breaking the symmetry somehow by using random
    // perturbation of a weak prior on exposures to try to improve mixing.
    // Another idea is to allow use of (e.g.) cosmic signatures as priors
    // on the signatures.
    for (s in 1:S) {
        exposures_prior[s] = exposures_prior_val;
    }
    
    for (c in 1:C) {
        signatures_prior[c] = 0.5;
    }
}
parameters {
    simplex[S] exposures[G];
    simplex[C] signatures[S];
}
transformed parameters {
    matrix<lower=0>[G, C] probs;
    {
        matrix[G, S] exposures_mat;
        matrix[S, C] signatures_mat;
    
        for (g in 1:G) {
            for (s in 1:S) {
                exposures_mat[g, s] = exposures[g, s];
            }
        }
        
        for (s in 1:S) {
            for (c in 1:C) {
                signatures_mat[s, c] = signatures[s, c];
            }
        }
        probs = exposures_mat * signatures_mat;
    }
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
        counts[g] ~ multinomial(to_vector(probs[g]));
    }
}
generated quantities {
    vector[G] log_lik;
    real bic;
    
    // Compute log_lik
    for (g in 1:G) {
        log_lik[g] = multinomial_lpmf(counts[g] | to_vector(probs[g]));
    }
    
    // Compute bic
    bic = 2 * sum(log_lik) - log(G) * (G*(S-1) + S*(C-1));
}
