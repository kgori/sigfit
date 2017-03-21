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
}
transformed parameters{
    matrix[G, C] probs;
    {
        matrix[G, S] exposures_mat;
        for (i in 1:G) {
            exposures_mat[i] = to_row_vector(exposures[i]);
        }
        probs = exposures_mat * signatures';
    }
}
model {
    for (i in 1:G) {
        exposures[i] ~ dirichlet(alpha);
        counts[i] ~ multinomial(to_vector(scale_to_sum_1(probs[i])));
    }
}
