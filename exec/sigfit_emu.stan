data {
    int N;          // Number of mutation channels [uses index j]
    int M;          // Number of tumour samples [uses index m]
    int n;          // Number of processes [uses index a]
    int X[N, M];    // Matrix of mutation counts
    int w[N, M];    // Matrix of opportunities
    vector[N] alpha;   // Prior on Dirichlet distribution for signatures
}
parameters {
    simplex[N] mu[n];           // Matrix of signatures, with simplex constraint
    matrix<lower=0>[n, M] x;    // Matrix of activities
}
transformed parameters {
    matrix<lower=0>[N, M] mu_x_w;    // elementwise product of mu_x and w
    {
        matrix[N, n] mu_mat;    // mu' recast as matrix
        matrix[N, M] mu_x;      // matrix product of mu and x
        for (j in 1:N) {
            for (a in 1:n) {
                mu_mat[j, a] = mu[a, j];
            }
        }
        mu_x = mu_mat * x;
        
        for (j in 1:N) {
            for (m in 1:M) {
                mu_x_w[j, m] = mu_x[j, m] * w[j, m];
            }
        }
    }
}
model {
    // Priors for x
    for (a in 1:n) {
        for (m in 1:M) {
            x[a, m] ~ cauchy(0, 2.5); // need something positive continuous
        }
    }
    
    // Priors for mu
    for (a in 1:n) {
        mu[a] ~ dirichlet(alpha);
    }
    
    // Likelihood
    for (j in 1:N) {
        X[j] ~ poisson(mu_x_w[j]);
    }
}
