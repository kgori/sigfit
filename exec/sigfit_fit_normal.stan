functions {
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
}
data {
    int G;
    int S;
    int C;
    matrix<lower=0>[G, C] counts;
    matrix<lower=0, upper=1>[S, C] signatures;
}
parameters {
    simplex[S] exposures[G];
}
transformed parameters {
    matrix[G, S] expos = array_to_matrix(exposures);
    matrix[G, C] R = expos * signatures;
}
model {
    for (g in 1:G) {
        exposures[g] ~ dirichlet(rep_vector(0.1, S));
        counts[g] ~ normal(R[g] * sum(counts[g]), 1);
    }
}
