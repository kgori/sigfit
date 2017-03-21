functions {
    vector enforcesum1(vector v) {
        return (v / sum(v));
    }
    row_vector enforcesum1row(row_vector v) {
        return (v / sum(v));
    }
}
data {
    int<lower=1> C; // number of categories
    int<lower=1> S; // number of signatures
    int<lower=1> G; // number of genomes
    matrix[C, S] signatures; // matrix of categories (rows) by signatures (columns)
    int counts[G, C]; // data = counts per category (columns) per genome sample (rows)
    vector<lower=0>[S] alpha; // prior on global probs
}
parameters {
    vector<lower=0>[G] kappa; // parameter to accommodate individual samples' departure from global trend
    simplex[S] globalExposures;
    simplex[S] sampleExposures[G];
}
transformed parameters {
    matrix<lower=0,upper=1>[G, C] probs;
    matrix[G, S] sampleExposuresMatrix;
    for (i in 1:G) {
        sampleExposuresMatrix[i] = to_row_vector(sampleExposures[i]);
    }
    probs = sampleExposuresMatrix * signatures';
    for (i in 1:G) {
        probs[i] = enforcesum1row(probs[i]);
    }
}
model {
    globalExposures ~ dirichlet(alpha);
    for (i in 1:G) {
        kappa[i] ~ cauchy(0, 2.5);
        sampleExposures[i] ~ dirichlet(kappa[G] * globalExposures);
        counts[i] ~ multinomial(to_vector(probs[i]));
    }
}
