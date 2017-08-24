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
