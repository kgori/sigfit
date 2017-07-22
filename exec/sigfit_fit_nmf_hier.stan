functions {
    #include "common_functions.stan"
}
data {
    int<lower=1> C; // number of categories
    int<lower=1> S; // number of signatures
    int<lower=1> G; // number of genomes
    matrix[C, S] signatures; // matrix of categories (rows) by signatures (columns)
    int counts[G, C]; // data = counts per category (columns) per genome sample (rows)
}
transformed data {
    vector[S] ones;  // Prior on shared weights
    for (i in 1:S) ones[i] = 1;
}
parameters {
    simplex[S] exposures[G];
    vector<lower=1,upper=1000>[S] alpha;
}
transformed parameters {
    vector<lower=0, upper=1>[C] probs[G];
    for (i in 1:G) {
        probs[i] = scale_to_sum_1(signatures * exposures[i]);
    }
}
model {
    alpha ~ uniform(1, 1000); // all exposures share this Dirichlet prior

    for (i in 1:G) {
        exposures[i] ~ dirichlet(alpha);
        counts[i] ~ multinomial(probs[i]);
    }
}
generated quantities {
    vector[G] log_lik;
    for (i in 1:G) {
        log_lik[i] = multinomial_lpmf(counts[i] | probs[i]);
    }
}
