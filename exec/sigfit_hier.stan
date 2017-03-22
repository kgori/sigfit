functions {
    vector scale_to_sum_1_v(vector v) {
        return (v / sum(v));
    }
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
    real<lower=0> shared_scale;
    simplex[S] shared_weights;
}
transformed parameters {
    vector[S] alpha = shared_scale * shared_weights;
}
model {
    vector[C] probs;

    shared_scale ~ cauchy(0, 2.5);
    shared_weights ~ dirichlet(ones);
    
    for (i in 1:G) {
        exposures[i] ~ dirichlet(alpha);
        probs = scale_to_sum_1_v(signatures * exposures[i]);
        counts[i] ~ multinomial(probs);
    }
}
generated quantities {
    vector[S] shared_exposures = dirichlet_rng(alpha);
}
