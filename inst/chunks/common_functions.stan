vector scale_to_sum_1(vector v) {
    return (v / sum(v));
}

row_vector scale_row_to_sum_1(row_vector r) {
    return (r / sum(r));
}
