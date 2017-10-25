functions {
    vector scale_to_sum_1(vector v) {
        return (v / sum(v));
    }
    
    row_vector scale_row_to_sum_1(row_vector r) {
        return (r / sum(r));
    }
    
    /**
       * Copy an array of equal-length vectors (or simplexes)
       * into a matrix
       *
       * @param x An array of vectors
       * @return A matrix copy of x
       */
    matrix array_to_matrix(vector[] x) {
        // Assume x doesn't have 0 rows or columns
        matrix[size(x), rows(x[1])] y;
        for (m in 1:size(x))
            y[m] = x[m]';
        return y;
    }
    //#include "common_functions.stan"
}
data {
    int<lower=1> C;            // number of mutation categories
    int<lower=1> S;            // number of mutational signatures
    int<lower=1> G;            // number of genomes
    matrix[S, C] signatures;   // matrix of signatures (columns) to be fitted
    int<lower=0> counts[G, C]; // data = counts per category (columns) per genome sample (rows)
    vector<lower=0>[S] alpha;  // prior on exposures (i.e. mixing proportions of signatures)
    matrix[G, C] opps;         // matrix of opportunities
}
parameters {
    simplex[S] exposures[G];
    real<lower=0> multiplier[G];
}
transformed parameters {
    matrix<lower=0>[G, S] exposures_raw;
    matrix[G, C] lambda;  // Poisson parameters
    
    for (g in 1:G) {
        exposures_raw[g] = exposures[g]' * multiplier[g];
    }

    // Poisson parameters
    lambda = exposures_raw * signatures .* opps;
}
model {
    for (i in 1:G) {
        // Priors
        exposures[i] ~ dirichlet(alpha);
        multiplier ~ cauchy(0, 1);
        
        // Likelihood
        counts[i] ~ poisson(lambda[i]);
    }
}
