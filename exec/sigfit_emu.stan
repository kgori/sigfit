data {
    int N;          // Number of mutation channels [uses index j]
    int M;          // Number of tumour samples [uses index m]
    int n;          // Number of processes [uses index a]
    int X[M, N];    // Matrix of mutation counts
    matrix[M, N] w;    // Matrix of opportunities
}
transformed data {
    vector[N] alpha;
    for (j in 1:N) {
        alpha[j] = 1;
    }
}
parameters {
    simplex[N] mu[n];           // Matrix of signatures, with simplex constraint
    matrix<lower=0>[M, n] x;
}
model {
    // Precomputation
    matrix[M, N] mu_x_w;    // elementwise product of mu_x and w
    matrix[n, N] mu_mat;    // mu recast as matrix
    
    for (a in 1:n) {
        for (j in 1:N) {
            mu_mat[a, j] = mu[a, j];
        }
    }
    mu_x_w = x * mu_mat .* w;

    // Priors for x
    for (m in 1:M) {
        x[m] ~ cauchy(0, 2.5); // need something positive continuous
    }
    
    // Priors for mu
    for (a in 1:n) {
        mu[a] ~ dirichlet(alpha);
    }
    
    // Likelihood
    for (m in 1:M) {
        X[m] ~ poisson(mu_x_w[m]);
    }
}
