functions {
#include /include/common_functions.stan
}

data {
    int<lower=1, upper=4> family;   // model: 1=multinomial, 2=poisson, 3=negbin, 4=normal
    int<lower=1> C;                 // number of mutation categories
    int<lower=1> S;                 // number of fixed signatures
    int<lower=1> G;                 // number of genomes
    int<lower=1> N;                 // number of extra signatures
    matrix[S, C] fixed_sigs;        // matrix of signatures (rows) by categories (columns)
    int counts_int[G, C];           // observed mutation counts (discrete case)
    real counts_real[G, C];         // observed mutation counts (continuous case)
    matrix[G, C] opportunities;     // mutational opportunities (genome per row)
    vector<lower=0>[S+N] kappa[G];  // prior on exposures (mixing proportions)
    vector<lower=0>[C] alpha[N];    // prior for extra signatures
    int<lower=0, upper=1> dpp;      // use DPP prior: 0=no, 1=yes
    real<lower=0> concentration;    // concentration hyperparameter of DPP
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
    simplex[C] extra_sigs[N];          // additional signatures to extract
    real<lower=0> multiplier[G_mult];  // exposure multipliers
    vector<lower=0>[G_sigma] sigma;    // standard deviations (normal/t model)
    vector<lower=0>[C_phi] phi;        // unscaled overdispersions (neg bin model)

    // Exposures parameter IF NOT using Dirichlet Process Prior
    simplex[T] exposures_all[G_notdpp];           // signature exposures (genome per row), for all signatures

    // Exposures priors IF using Dirichlet Process Prior
    simplex[S] exposures_fixed[G_dpp];                               // Exposures within the fixed set
    vector<lower=0, upper=1>[N + 1] exposures_sticklengths[G_dpp];   // signature exposures (genome per row)
    real<lower=0> dp_alpha[G_dpp];                                   // dirichlet process prior
}

transformed parameters {
    simplex[T] exposures[G];
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
        multiplier ~ cauchy(0, 2.5);

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
