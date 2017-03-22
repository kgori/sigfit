functions {
    vector scale_to_sum_1(vector v) {
        return (v / sum(v));
    }
}
data {
    int<lower=1> C; // number of mutation categories
    int<lower=1> S; // number of mutational signatures
    int<lower=1> G; // number of genomes
    matrix[C, S] signatures; // matrix: each row is probability vec. of generating each cat. of mutation
    int counts[G, C]; // data = counts per category (columns) per genome sample (rows)
    vector<lower=0>[S] alpha; // prior on exposures (i.e. mixing proportions of signatures)
}
parameters {
    simplex[S] exposures[G];
}
model {
    vector[C] probs;
    for (i in 1:G) {
        exposures[i] ~ dirichlet(alpha);
        probs = scale_to_sum_1(signatures * exposures[i]);
        counts[i] ~ multinomial(probs);
    }
}
