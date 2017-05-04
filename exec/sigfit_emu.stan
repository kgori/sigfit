data {
    int N;          // Number of mutation channels [uses index j]
    int M;          // Number of tumour samples [uses index m]
    int n;          // Number of signatures [uses index a]
    int counts[M, N];    // Matrix of mutation counts
    matrix[M, N] opps; // Matrix of opportunities
}
transformed data {
    vector[N] alpha;
    for (j in 1:N) {
        alpha[j] = 1;
    }
}
parameters {
    simplex[N] signatures[n];           // Matrix of signatures, with simplex constraint
    matrix<lower=0>[M, n] exposures;
}
transformed parameters {
    // Precomputation
    matrix[M, N] sig_expo_opps;    // elementwise product of sig_expo and opps
    matrix[n, N] sig_mat;    // sig recast as matrix
    
    for (a in 1:n) {
        for (j in 1:N) {
            sig_mat[a, j] = signatures[a, j];
        }
    }
    sig_expo_opps = exposures * sig_mat .* opps;
}
model {
    // Priors for x
    for (m in 1:M) {
        exposures[m] ~ cauchy(0, 2.5); // need something positive continuous
    }
    
    // Priors for mu
    for (a in 1:n) {
        signatures[a] ~ dirichlet(alpha);
    }
    
    // Likelihood
    for (m in 1:M) {
        counts[m] ~ poisson(sig_expo_opps[m]);
    }
}
generated quantities {
    matrix[M, N] log_lik;
    for (m in 1:M) {
        for (j in 1:N) {
            log_lik[m, j] = poisson_lpmf(counts[m, j] | sig_expo_opps[m, j]);
        }
    }
}
