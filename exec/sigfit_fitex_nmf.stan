data {
    int<lower=1> C;   // number of categories
    int<lower=1> S;   // number of fixed signatures
    int<lower=1> G;   // number of genomes
    int<lower=1> N;   // number of extra signatures
    matrix[S, C] fixed_sigs;  // matrix of signatures (rows) by categories (columns)
    int counts[G, C];         // data = counts per genome (rows) in each category (columns)
}
transformed data {
    int T = S + N;              // total number of signatures, including extra signatures
    vector[C] signature_prior;  // Jeffreys prior for extra_signature
    vector[T] exposures_prior;  // Jeffreys prior for exposures
    for (c in 1:C) {
        signature_prior[c] = 0.5;
    }
    for (s in 1:T) {
        exposures_prior[s] = 0.5;
    }
}
parameters {
    simplex[C] extra_sigs[N];  // additional signatures to extract
    simplex[T] exposures[G];   // includes exposure for extra_signature
}
transformed parameters {
    matrix[N, C] extra_sigs_mat;
    matrix[G, T] exposures_mat;
    matrix[T, C] signatures;
    matrix<lower=0>[G, C] probabilities;
    for (i in 1:G) {
        for (j in 1:T) {
            exposures_mat[i, j] = exposures[i, j];
        }
    }
    for (i in 1:N) {
        for (j in 1:C) {
            extra_sigs_mat[i, j] = extra_sigs[i, j];
        }
    }
    signatures = append_row(fixed_sigs, extra_sigs_mat);
    probabilities = exposures_mat * signatures;
}
model {
    for (n in 1:N) {
        extra_sigs[n] ~ dirichlet(signature_prior);
    }
    for (g in 1:G) {
        exposures[g] ~ dirichlet(exposures_prior);
        counts[g] ~ multinomial(to_vector(probabilities[g]));
    }
}
generated quantities {
    vector[G] log_lik;
    for (g in 1:G) {
        log_lik[g] = multinomial_lpmf(counts[g] | to_vector(probabilities[g]));
    }
}
