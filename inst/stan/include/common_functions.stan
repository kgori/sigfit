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

/**
    * https://ecosang.github.io/blog/study/dirichlet-process-with-stan/
    * Generates a simplex using a stick-breaking method.
    * @param v A vector of proportions (values between 0 and 1)
    * @return A vector that sums to 1
    */
vector stick_breaking(vector v) {
    int C = num_elements(v);
    vector[C] pi = rep_vector(0, C);
    pi[1] = v[1];
    for (j in 2:C) {
        pi[j]= v[j]*(1-v[j-1])*pi[j-1]/v[j-1];
    }
    return scale_to_sum_1(pi);
}

/**
    * Helper function for fit-extract dpp model.
    * Returns a simplex of size T (T=N+S) from two exposures vectors: sticks (size N+1)
    * and fixed (size S). Here, S is the number of fixed signatures and N is the number
    * of additional extracted signatures. `sticks` is sampled from the DPP and `fixed`
    * is sampled from a straight Dirichlet.
    *
    * Why is sticks length N+1?
    * The first element of `sticks` is the weight the model gives to all fixed
    * signatures put together; the remaining elements are the weights of the additional
    * signatures. The first weight is split according to the weights in `fixed`. This
    * results in a vector of weights for all signatures combined.
    */
vector build_exposures_from_sticks_and_fixed(vector sticks, vector fixed) {
    int N = num_elements(sticks) - 1; // N = number of extra signatures
    int S = num_elements(fixed);      // S = number of fixed signatures
    int T = N + S;                    // T = total number of signatures
    vector[N+1] dpp_weights;
    vector[S] fixed_weights;
    vector[T] exposures;
    dpp_weights = stick_breaking(sticks);
    fixed_weights = dpp_weights[1] * fixed; // subdivide DPP weight 1 by Dirichlet weights in `fixed`

    // Put weights in place in exposures
    exposures[1:S] = fixed_weights;
    exposures[(S+1):T] = dpp_weights[2:(N+1)];

    // Should sum to 1 already, but let's make sure
    return scale_to_sum_1(exposures);
}
