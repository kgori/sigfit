data {
    int<lower=0> T;
    int<lower=0> I;
    int<lower=0> K;
    matrix<lower=0.0>[T,I] counts;
    real<lower=0> sigma[I];
}
transformed data {
    real<lower=0> e_bar;
    real<lower=0> e_sigma;
    vector[T] temp;
    for (t in 1:T) 
        temp[t] = log(sum(counts[t]));
    e_bar = mean(temp);
    e_sigma = sd(temp);
}
parameters {
    matrix<lower=0>[T,K] exposures;
    simplex[I] signatures[K];
}
transformed parameters {
    matrix<lower=0>[T,I] mu;
    {
        matrix[K,I] signatures_mat;
        for (k in 1:K) {
            signatures_mat[k] = to_row_vector(signatures[k]);
        }
        mu = exposures * signatures_mat;
    }
}
model {
    // Prior
    for (t in 1:T)
        exposures[t] ~ lognormal(e_bar, e_sigma);
    
    // Likelihood
    for (t in 1:T) {
        for (i in 1:I) {
            counts[t, i] ~ normal(mu[t, i], sigma[i]) T[0, ];
        }
    }
}
generated quantities {
    vector[T*I] log_lik;
    for (t in 1:T) {
        for (i in 1:I) {
            // see stan manual p79 for pdf of truncated normal
            if (counts[t, i] < 0)
                log_lik[(t - 1) * I + i] = negative_infinity();
            else
                log_lik[(t - 1) * I + i] = normal_lpdf(counts[t, i] | mu[t, i], sigma[i]) - normal_lccdf(counts[t, i] | mu[t, i], sigma[i]);
        }
    }
}
