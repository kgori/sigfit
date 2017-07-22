data {
    int N;  // Number of mutation channels [uses index j]
    int M;  // Number of tumour samples [uses index m]
    int n;  // Number of processes [uses index a]
    int X[N, M];       // Matrix of mutation counts
    int w[N, M];       // Matrix of opportunities
    simplex[N] mu[n];  // Matrix of signatures
}
transformed data {
    # Signatures and opportunities recast as matrices
    matrix[N, n] mu_mat;
    matrix[N, M] w_mat;
    
    for (a in 1:n) {
        for (j in 1:N) {
            mu_mat[j,a] = mu[a,j];
        }
    }
    
    for (j in 1:N) {
        for (m in 1:M) {
            w_mat[j,m] = w[j,m];
        }
    }
}
parameters {
    matrix<lower=0>[n, M] exposures;    // Matrix of exposures (activities)
}
model {
    matrix[N, M] mu_expo_w;
    
    // Priors for exposures
    for (a in 1:n) {
        for (m in 1:M) {
            exposures[a, m] ~ cauchy(0, 2.5); // need something positive continuous
        }
    }
    
    mu_expo_w = mu_mat * exposures .* w_mat;
    
    // Likelihood
    for (j in 1:N) {
        X[j] ~ poisson(mu_expo_w[j]);
    }
}
