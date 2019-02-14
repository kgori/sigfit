functions {
//    #include "common_functions.stan"
    vector scale_to_sum_1(vector v) {
        return (v / sum(v));
    }
    row_vector scale_row_to_sum_1(row_vector r) {
        return (r / sum(r));
    }
    matrix array_to_matrix(vector[] x) {
        // Assume x doesn't have 0 rows or columns
        matrix[size(x), rows(x[1])] y;
        for (m in 1:size(x))
            y[m] = x[m]';
        return y;
    }
}
data {
    int<lower=1> C;               // number of mutation categories
    int<lower=1> S;               // number of mutational signatures
    int<lower=1> G;               // number of genomes
    vector[C] counts[G];          // observed mutation counts (row per genome)
    matrix[S, C] signatures;      // signatures to fit (row per signature)
}
parameters {
    vector<lower=0, upper=pi()/2>[S] activities_unif[G];
    vector<lower=0, upper=pi()/2>[C] nu_unif[G];
    vector<lower=0, upper=pi()/2>[G] sigma_unif;
}
transformed parameters {
    vector<lower=0>[S] activities[G];   // signature activities (row per genome)
    vector<lower=0>[C] nu[G];           // degrees of freedom
    vector<lower=0>[G] sigma;           // scales
    matrix[G, C] expected_counts;
    activities = tan(activities_unif);  // Reparameterisation of Half-Cauchy(0, 1)
    nu = tan(nu_unif);                  //  from U(0, pi/2)
    sigma = tan(sigma_unif);
    expected_counts = array_to_matrix(activities) * signatures;
}
model {
    for (g in 1:G) {
        counts[g] ~ student_t(nu[g], expected_counts[g], sigma[g]);
    }
}
generated quantities {
    vector[S] exposures[G];
    for (g in 1:G) {
        exposures[g] = scale_to_sum_1(activities[g]);
    }
}
