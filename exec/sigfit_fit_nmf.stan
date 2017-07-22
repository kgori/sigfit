functions {
    #include "common_functions.stan"
}
data {
    int<lower=1> C; // number of mutation categories
    int<lower=1> S; // number of mutational signatures
    int<lower=1> G; // number of genomes
    matrix[C, S] signatures; // matrix of signatures (columns) to be fitted
    int counts[G, C]; // data = counts per category (columns) per genome sample (rows)
    vector<lower=0>[S] alpha; // prior on exposures (i.e. mixing proportions of signatures)
}
parameters {
    simplex[S] exposures[G];
}
transformed parameters {
    vector<lower=0, upper=1>[C] probs[G];
    for (i in 1:G) {
        probs[i] = scale_to_sum_1(signatures * exposures[i]);
    }
}
model {
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
