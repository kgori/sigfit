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
    vector[C] pi;
    pi[1] = v[1];
    // stick-break process based on The BUGS book Chapter 11 (p.294)
    for(j in 2:(C-1)){
        pi[j] = v[j] * (1 - v[j-1]) * pi[j-1] / v[j-1];
    }
    pi[C] = 1 - sum(pi[1:(C-1)]); // to make a simplex
    return pi;
}
