data {
    int<lower=1> G;
    int<lower=1> C;
    int<lower=1> S;
    int<lower=0> counts[G, C];
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
    matrix<lower=0>[G, C] probabilities;
    {
        matrix[G, S] exposuresMat;
        matrix[S, C] signaturesMat;
    
        for (i in 1:G) {
            for (j in 1:S) {
                exposuresMat[i, j] = exposures[i, j];
            }
        }
        
        for (i in 1:S) {
            for (j in 1:C) {
                signaturesMat[i, j] = signatures[i, j];
            }
        }
        probabilities = exposuresMat * signaturesMat;
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
        counts[g] ~ multinomial(to_vector(probabilities[g]));
    }
}
generated quantities {
    vector[G] log_lik;
    for (g in 1:G) {
        log_lik[g] = multinomial_lpmf(counts[g] | to_vector(probabilities[g]));
    }
}
