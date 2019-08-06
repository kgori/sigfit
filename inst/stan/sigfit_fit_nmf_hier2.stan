functions {
#include /include/common_functions.stan
}
data {
    int<lower=1> C; // number of categories
    int<lower=1> S; // number of signatures
    int<lower=1> G; // number of genomes
    int<lower=1> K; // number of groups
    int<lower=1,upper=K> membership[G]; // group membership of each sample
    matrix[S, C] signatures; // matrix of categories (rows) by signatures (columns)
    int counts[G, C]; // data = counts per category (columns) per genome sample (rows)
    real<lower=0> gamma_alpha; // value of gamma prior alpha parameter (4 is good)
    real<lower=0> gamma_beta; // value of gamma prior beta parameter (0.4 is good)
}
parameters {
    simplex[S] exposures[G];
    vector<lower=0>[S] group_alphas[K];
}
transformed parameters {
    matrix<lower=0, upper=1>[G, C] probs = array_to_matrix(exposures) * signatures;
}
model {
    for (k in 1:K) {
        group_alphas[k] ~ gamma(gamma_alpha, gamma_beta);
    }
    for (g in 1:G) {
        exposures[g] ~ dirichlet(group_alphas[membership[g]]);
        counts[g] ~ multinomial(probs[g]');
    }
}
generated quantities {
    simplex[S] group_exposures[K];
    for (k in 1:K) {
        group_exposures[k] = scale_to_sum_1(group_alphas[k]);
    }
}
