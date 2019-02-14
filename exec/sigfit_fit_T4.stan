functions {
//    #include "common_functions.stan"
    vector scale_to_sum_1(vector v) {
        return (v / sum(v));
    }
    row_vector scale_row_to_sum_1(row_vector r) {
        return (r / sum(r));
    }
}
data {
    int<lower=1> C;               // number of mutation categories
    int<lower=1> S;               // number of mutational signatures
    int<lower=1> G;               // number of genomes
    matrix[G, C] counts;          // observed mutation counts (row per genome)
    matrix[S, C] signatures;      // signatures to fit (row per signature)
}
parameters {
    matrix<lower=0>[G, S] activities;  // signature activities (row per genome)
    vector<lower=1>[C] nu;             // degrees of freedom
    vector<lower=0>[G] sigma;          // standard deviations
}
transformed parameters {
    matrix[G, C] expected_counts = activities * signatures;
}
model {
    sigma ~ cauchy(0, 1);
    nu ~ cauchy(0, 1);
    for (g in 1:G) {
        activities[g] ~ cauchy(0, 1);
        counts[g] ~ student_t(nu, expected_counts[g], sigma[g]);
    }
}
generated quantities {
    matrix[G, S] exposures;
    for (g in 1:G) {
        exposures[g] = scale_row_to_sum_1(activities[g]);
    }
}
