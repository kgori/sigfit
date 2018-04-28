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

//    #include "common_functions.stan"
}
data {
    int<lower=1> C; // number of categories
    int<lower=1> S; // number of signatures
    int<lower=1> G; // number of genomes
    int<lower=1> K; // number of groups
    int<lower=1,upper=K> membership[G]; // group membership of each sample
    matrix[S, C] signatures; // matrix of categories (rows) by signatures (columns)
    int counts[G, C]; // data = counts per category (columns) per genome sample (rows)
    vector<lower=0>[S] alpha; // prior on group exposures
}
parameters {
    vector[G] genome_weights_raw;
}
transformed parameters {
    simplex[S] group_exposures[K];
    matrix<lower=0,upper=1>[K, C] probs;
    probs = array_to_matrix(group_exposures) * signatures;
}
model {
    real group_sums[K];
    genome_weights_raw ~ gamma(1, 1);
    for (k in 1:K) group_sums[k] = 0;
}
