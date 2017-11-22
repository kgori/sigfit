functions {
    #include "common_functions.stan"
}
data {
    int<lower=1> C; // number of categories
    int<lower=1> S; // number of signatures
    int<lower=1> G; // number of genomes
    int<lower=1> K; // number of groups
    int<lower=1,upper=K> membership[G]; // group membership of each sample
    matrix[S, C] signatures; // matrix of categories (rows) by signatures (columns)
    int counts[G, C]; // data = counts per category (columns) per genome sample (rows)
    vector<lower=0>[S] alpha; // prior on group exposures
}
parameters {
    vector<lower=0>[K] concentration; // parameter to accommodate individual samples' departure from global trend
    simplex[S] group_exposures[K];
    simplex[S] exposures[G];
}
transformed parameters {
    matrix<lower=0,upper=1>[G, C] probs;
    probs = array_to_matrix(exposures) * signatures;
}
model {
    for (k in 1:K) {
        group_exposures[k] ~ dirichlet(alpha);
        concentration[k] ~ cauchy(0, 2.5);
    }
    for (g in 1:G) {
        exposures[g] ~ dirichlet(concentration[membership[g]] * 
            group_exposures[membership[g]]);
        counts[g] ~ multinomial(probs[g]');
    }
}
