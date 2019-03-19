functions {
    #include "common_functions.stan"
}

data {
    int<lower=1, upper=3> family;  // model: 1=multinomial, 2=poisson, 3=normal
    int<lower=0, upper=1> robust;  // robust model: 0=no, 1=yes (neg binomial or t)
    int<lower=1> C;                // number of mutation categories
    int<lower=1> S;                // number of mutational signatures
    int<lower=1> G;                // number of genomes
    int counts_int[G, C];          // observed mutation counts (discrete case)
    real counts_real[G, C];        // observed mutation counts (continuous case)
    matrix[S, C] signatures;       // signatures to fit (signature per row)
    matrix[G, C] opportunities;    // mutational opportunities (genome per row)
    vector<lower=0>[S] kappa;      // prior on exposures (mixing proportions)
}

transformed data {
    // Dynamic dimensions for model-specific parameters:
    // unused parameters have zero length
    int C_phi = ((family == 2) && (robust == 1)) ? C : 0;
    int C_nu = ((family == 3) && (robust == 1)) ? C : 0;
    int G_sigma = (family == 3) ? G : 0;
    int G_mult = (family != 1) ? G : 0;
}

parameters {
    simplex[S] exposures[G];           // signature exposures (genome per row)
    real<lower=0> multiplier[G_mult];  // exposure multipliers
    vector<lower=0>[G_sigma] sigma;    // standard deviations (normal/t model)
    vector<lower=1>[C_nu] nu;          // degrees of freedom (t model)
    vector<lower=0>[C_phi] phi_raw;    // unscaled overdispersions (neg bin model)
}

transformed parameters {
    matrix<lower=0>[G, S] activities;  // scaled exposures (# mutations)
    matrix[G, C] expected_counts;
    vector<lower=0>[C_phi] phi;
    // Scale overdispersion (neg binomial model)
    if ((family == 2) && (robust == 1)) {
        for (c in 1:C) {
            phi[c] = 1 / phi_raw[c] ^ 2;
        }
    }
    // Scale exposures into activities
    if (family == 1) {
        // Multinomial model uses unscaled exposures
        activities = array_to_matrix(exposures);
    }
    else if (family == 2) {
        for (g in 1:G) {
            activities[g] = exposures[g]' * sum(counts_int[g]) * multiplier[g];
        }
    }
    else {
        for (g in 1:G) {
            activities[g] = exposures[g]' * sum(counts_real[g]) * multiplier[g];
        }
    }
    // Calculate expected counts (or probabilities)
    expected_counts = activities * signatures .* opportunities;
}

model {
    // Exposure priors (all models)
    for (g in 1:G) {
        exposures[g] ~ dirichlet(kappa);
    }
    
    // Multinomial ('NMF') model
    if (family == 1) {
        for (g in 1:G) {
            counts_int[g] ~ multinomial(expected_counts[g]');
        }
    }
    
    else {
        multiplier ~ cauchy(0, 2.5);
        
        // Poisson model family
        if (family == 2) {
            
            // Poisson ('EMu') model
            if (robust == 0) {
                for (g in 1:G) {
                    counts_int[g] ~ poisson(expected_counts[g]);
                }
            }
            
            // Negative binomial model
            else {
                phi_raw ~ normal(0, 1);
                for (g in 1:G) {
                    counts_int[g] ~ neg_binomial_2(expected_counts[g], phi);
                }
            }
        }
    
        // Normal model family
        else if (family == 3) {
            sigma ~ cauchy(0, 2.5);
            
            // Normal model
            if (robust == 0) {
                for (g in 1:G) {
                    counts_real[g] ~ normal(expected_counts[g], sigma[g]);
                }
            }
            
            // t model
            else {
                nu ~ gamma(2, 0.1);
                for (g in 1:G) {
                    counts_real[g] ~ student_t(nu, expected_counts[g], sigma[g]);
                }
            }
        }
    }
}
