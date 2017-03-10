data {
    int<lower=1> C; // number of categories
    int<lower=1> S; // number of signatures
    matrix[C, S] signatures; // matrix of categories (rows) by signatures (columns)
    int counts[C]; // data = counts per category
    vector<lower=0>[S] alpha; // prior on probs
}
parameters {
    simplex[S] exposures;
}
transformed parameters {
    simplex[C] probs;
    probs = signatures * exposures;
    probs = probs / sum(probs);
}
model {
    exposures ~ dirichlet(alpha);
    counts ~ multinomial(probs);
}

