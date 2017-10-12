functions{
    #include "common_functions.stan"
}
data {
    int<lower=1> C;             // number of mutation categories
    int<lower=1> S;             // number of mutational signatures
    int<lower=1> G;             // number of genomes
    int<lower=0> counts[G, C];  // matrix of counts per category (columns) per genome sample (rows)
    matrix[S, C] alpha;         // prior for signatures
}
transformed data {
    // Not put much thought into this prior - just used the Jeffreys prior.
    // One thing to try would be breaking the symmetry somehow by using random
    // perturbation of a weak prior on exposures to try to improve mixing.
    vector[S] kappa = rep_vector(0.5, S);
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
    // Currently nothing is being done about it in the model,
    // and probably avoiding running multiple chains will be fine.
    
    for (s in 1:S) {
        // Priors for signatures
        signatures[s] ~ dirichlet(alpha[s]');
    }
    
    for (g in 1:G) {
        // Priors for exposures (Jeffreys)
        exposures[g] ~ dirichlet(kappa);
        
        // Likelihood
        counts[g] ~ multinomial(probs[g]');
    }
}
generated quantities {
    vector[G] log_lik;
    real bic;
    
    // Compute log likelihood
    for (g in 1:G) {
        log_lik[g] = multinomial_lpmf(counts[g] | probs[g]');
    }
    
    // Compute BIC with G*(S-1) + S*(C-1) parameters
    bic = 2 * sum(log_lik) - log(G) * (G*(S-1) + S*(C-1));
}
