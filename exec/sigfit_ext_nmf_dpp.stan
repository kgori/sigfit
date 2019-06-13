functions{
    #include "common_functions.stan"
}
data {
    int<lower=1> C;             // number of mutation categories
    int<lower=1> S;             // number of mutational signatures
    int<lower=1> G;             // number of genomes
    int<lower=0> counts[G, C];  // matrix of counts per category (columns) per genome sample (rows)
    matrix[S, C] alpha;         // prior for signatures
    real<lower=0> kappa;        // prior for exposures
}
parameters {
    real<lower=0> dp_alpha[G]; // dirichlet process concentration hyperparameter
    vector<lower=0, upper=1>[S] exposures_sticklengths[G];
    simplex[C] signatures[S];
}
transformed parameters {
    matrix<lower=0, upper=1>[G, C] probs;
    simplex[S] exposures[G];
    for (g in 1:G) {
        exposures[g] = stick_breaking(exposures_sticklengths[g]);
    }
    // array_to_matrix is defined in common_functions.stan and is not in base Stan
    probs = array_to_matrix(exposures) * array_to_matrix(signatures);
}
model {
    dp_alpha ~ gamma(kappa, 1);

    for (s in 1:S) {
        // Priors for signatures
        signatures[s] ~ dirichlet(alpha[s]');
    }

    for (g in 1:G) {
        // Priors for exposures uniform dirichlet
        exposures_sticklengths[g] ~ beta(1, dp_alpha[g]);

        // Likelihood
        counts[g] ~ multinomial(probs[g]');
    }
}
