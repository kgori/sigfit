functions {
#include /include/common_functions.stan
}

data {
    int<lower=1, upper=4> family;         // model: 1=multinomial, 2=poisson, 3=negbin, 4=normal
    int<lower=1> C;                       // number of mutation categories
    int<lower=1> S;                       // number of fixed signatures
    int<lower=1> G;                       // number of genomes
    int<lower=1> N;                       // number of extra signatures
    matrix[S, C] fixed_sigs;              // matrix of signatures (rows) by categories (columns)
    array[G, C] int counts_int;           // observed mutation counts (discrete case)
    array[G, C] real counts_real;         // observed mutation counts (continuous case)
    matrix[G, C] opportunities;           // mutational opportunities (genome per row)
    array[G] vector<lower=0>[S+N] kappa;  // prior on exposures (mixing proportions)
    array[N] vector<lower=0>[C] alpha;    // prior for extra signatures
    int<lower=0, upper=1> dpp;            // use DPP prior: 0=no, 1=yes
    real<lower=0> concentration;          // concentration hyperparameter of DPP
}

transformed data {
    // Dynamic dimensions for model-specific parameters:
    // unused parameters have zero length
    int G_mult = (family != 1) ? G : 0;
    int C_phi = (family == 3) ? C : 0;
    int G_sigma = (family == 4) ? G : 0;
    int G_dpp = (dpp == 1) ? G : 0;
    int G_notdpp = (dpp == 0) ? G : 0;
    int T = S + N;  // total number of signatures
}

parameters {
    array[N] simplex[C] extra_sigs;    // additional signatures to extract
    array[G_mult] real<lower=0> mult;  // exposure multipliers
    vector<lower=0>[G_sigma] sigma;    // standard deviations (normal/t model)
    vector<lower=0>[C_phi] phi;        // unscaled overdispersions (neg bin model)

    // Exposures parameter IF NOT using Dirichlet Process Prior
    array[G_notdpp] simplex[T] exposures_all;  // signature exposures (genome per row), for all signatures

    // Exposures priors IF using Dirichlet Process Prior
    array[G_dpp] simplex[S] exposures_fixed;                              // Exposures within the fixed set
    array[G_dpp] vector<lower=0, upper=1>[N + 1] exposures_sticklengths;  // signature exposures (genome per row)
    array[G_dpp] real<lower=0> dp_alpha;                                  // dirichlet process prior
}

transformed parameters {
    array[G] simplex[T] exposures;
    matrix<lower=0>[G, T] activities;  // scaled exposures (# mutations)
    matrix[G, C] expected_counts;
    matrix[T, C] signatures = append_row(fixed_sigs, array_to_matrix(extra_sigs));

    // Construct the exposures
    if (dpp == 1) {
        for (g in 1:G) {
            exposures[g] = build_exposures_from_sticks_and_fixed(exposures_sticklengths[g], exposures_fixed[g]);
        }
    }
    else {
        for (g in 1:G) {
            exposures[g] = exposures_all[g];
        }
    }

    // Scale exposures into activities
    if (family == 1) {
        // Multinomial model uses unscaled exposures
        activities = array_to_matrix(exposures);
    }
    else if ((family == 2) || (family == 3)) {
        for (g in 1:G) {
            activities[g] = exposures[g]' * sum(counts_int[g]) * mult[g];
        }
    }
    else {
        for (g in 1:G) {
            activities[g] = exposures[g]' * sum(counts_real[g]) * mult[g];
        }
    }
    // Calculate expected counts (or probabilities)
    expected_counts = activities * signatures .* opportunities;

    if (family == 1) { // Multinomial requires a simplex
        for (g in 1:G) {
            expected_counts[g] = scale_row_to_sum_1(expected_counts[g]);
        }
    }
}

model {
    if (dpp == 1) {
        dp_alpha ~ gamma(concentration, 1);
    }

    // Exposure priors (all models)
    for (g in 1:G) {
        if (dpp == 1) {
            exposures_fixed[g] ~ dirichlet(rep_vector(1, S));
            exposures_sticklengths[g] ~ beta(1, dp_alpha[g]);
        }
        else {
            exposures_all[g] ~ dirichlet(kappa[g]);
        }
    }

    for (n in 1:N) {
        // Priors for signatures
        extra_sigs[n] ~ dirichlet(alpha[n]);
    }

    // Multinomial ('NMF') model
    if (family == 1) {
        for (g in 1:G) {
            counts_int[g] ~ multinomial(expected_counts[g]');
        }
    }

    else {
        mult ~ cauchy(0, 2.5);

        // Poisson model
        if (family == 2) {
            for (g in 1:G) {
                counts_int[g] ~ poisson(expected_counts[g]);
            }
        }

        // Negative binomial model
        else if (family == 3) {
            phi ~ cauchy(0, 2.5);
            for (g in 1:G) {
                counts_int[g] ~ neg_binomial_2(expected_counts[g], phi);
            }
        }

        // Normal model
        else {
            sigma ~ cauchy(0, 2.5);
            for (g in 1:G) {
                counts_real[g] ~ normal(expected_counts[g], sigma[g]);
            }
        }
    }
}
