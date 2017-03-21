functions {
    row_vector scale_to_sum_1(row_vector v) {
        return (v / sum(v));
    }
}
data {
    int<lower=1> C; // number of categories
    int<lower=1> S; // number of signatures
    int<lower=1> G; // number of genomes
    matrix[C, S] signatures; // matrix of categories (rows) by signatures (columns)
    int counts[G, C]; // data = counts per category (columns) per genome sample (rows)
    vector<lower=0>[S] alpha; // prior on probs
}
parameters {
    simplex[S] exposures[G];
    vector<lower=0>[S] shared_prior;
}
transformed parameters{
    matrix<lower=0,upper=1>[G, C] probs;
    {
        matrix[G, S] exposures_mat;
        for (i in 1:G) {
            exposures_mat[i] = to_row_vector(exposures[i]);
        }
        probs = exposures_mat * signatures';
    }
}
model {
    for (i in 1:S) {
        shared_prior[S] ~ cauchy(0, 2.5) T[, 100];
    }
    for (i in 1:G) {
        exposures[i] ~ dirichlet(shared_prior);
        counts[i] ~ multinomial(to_vector(scale_to_sum_1(probs[i])));
    }
}
generated quantities {
    vector[S] global_exposures = dirichlet_rng(shared_prior);
}
