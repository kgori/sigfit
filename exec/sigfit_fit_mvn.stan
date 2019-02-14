// functions {
//     //#include "common_functions.stan"
//     vector scale_to_sum_1(vector r) {
//         return (r / sum(r));
//     }
// }
// data {
//     int<lower=1> C;            // number of mutation categories
//     int<lower=1> S;            // number of mutational signatures
//     matrix[S, C] signatures;   // matrix of signatures (rows) to be fitted
//     vector<lower=0>[C] counts; // ** G=1
//     vector<lower=0>[S] kappa;  // prior on exposures (i.e. mixing proportions of signatures)
// }
// parameters {
//     //vector<lower=0>[S] activities;
//     simplex[S] exposures;
//     real<lower=0> multiplier;
//     cholesky_factor_corr[C] L_Omega;
//     vector<lower=0>[C] sigma;       // standard deviations
//     vector[C] alpha;
// }
// transformed parameters {
//     vector<lower=0>[S] activities = exposures * multiplier;
//     //matrix[C, C] L_Sigma = diag_pre_multiply(sigma, L_Omega);  // Cholesky factor of covariance
// }
// model {
//     vector[C] mu = to_vector(activities' * signatures);  // means of MVN distribution
//     exposures ~ dirichlet(kappa);
//     multiplier ~ cauchy(0, 2.5);
//     sigma ~ cauchy(0, 2.5);
//     L_Omega ~ lkj_corr_cholesky(3);
//     alpha ~ normal(0, 1);  // implies counts ~ multi_normal_cholesky(expected_counts, L_Sigma)
//     mu = counts - sigma .* (L_Omega * alpha);  // counts = mu + diag_matrix(sigma)*L*alpha
// }
// // generated quantities {
// //     vector<lower=0, upper=1>[S] exposures = scale_to_sum_1(activities);
// // }



functions {
    //#include "common_functions.stan"
    matrix array_to_matrix(vector[] x) {
         // Assume x doesn't have 0 rows or columns
         matrix[size(x), rows(x[1])] y;
         for (m in 1:size(x))
             y[m] = x[m]';
         return y;
    }
}
data {
    int<lower=1> C;            // number of mutation categories
    int<lower=1> S;            // number of mutational signatures
    int<lower=1> G;            // number of genomes
    matrix[S, C] signatures;   // matrix of signatures (rows) to be fitted
    matrix[G, C] counts;       // matrix of counts per category (columns) per sample (rows)
    vector<lower=0>[S] kappa;  // prior on exposures (i.e. mixing proportions of signatures)
}
parameters {
    simplex[S] exposures[G];             // signature exposures (mixing proportions)
    real<lower=0> multiplier[G];
    cholesky_factor_corr[C] L_Omega[G];  // Cholesky factor of correlation (Stan manual, S9.15)
    vector<lower=0>[C] L_sigma[G];       // standard deviations
}
transformed parameters {
    matrix<lower=0>[G, S] activities;
    matrix<lower=0>[G, C] mu;          // means of MVN distribution
    matrix[C, C] L_Sigma[G];           // Cholesky factor of covariance
    for (g in 1:G) {
        activities[g] = exposures[g]' * multiplier[g];
        L_Sigma[g] = diag_pre_multiply(L_sigma[g], L_Omega[g]);
    }
    mu = activities * signatures;
}
model {
    // Priors
    multiplier ~ cauchy(0, 1);
    for (g in 1:G) {
        exposures[g] ~ dirichlet(kappa);
        L_Omega[g] ~ lkj_corr_cholesky(3);
        L_sigma[g] ~ cauchy(0, 2.5);

        // Likelihood
        counts[g] ~ multi_normal_cholesky(mu[g], L_Sigma[g]);
    }
}




// functions {
//     //#include "common_functions.stan"
//     matrix array_to_matrix(vector[] x) {
//          // Assume x doesn't have 0 rows or columns
//          matrix[size(x), rows(x[1])] y;
//          for (m in 1:size(x))
//              y[m] = x[m]';
//          return y;
//     }
// }
// data {
//     int<lower=1> C;            // number of mutation categories
//     int<lower=1> S;            // number of mutational signatures
//     int<lower=1> G;            // number of genomes
//     matrix[S, C] signatures;   // matrix of signatures (rows) to be fitted
//     matrix[G, C] counts;       // matrix of counts per category (columns) per sample (rows)
//     vector<lower=0>[S] kappa;  // prior on exposures (i.e. mixing proportions of signatures)
// }
// parameters {
//     simplex[S] exposures[G];             // signature exposures (mixing proportions)
//     real<lower=0> multiplier[G];
//     cholesky_factor_corr[C] L_Omega[G];  // Cholesky factor of correlation (Stan manual, S9.15)
//     matrix<lower=0>[G, C] L_sigma;       // standard deviations
// }
// transformed parameters {
//     matrix<lower=0>[G, S] activities;
//     matrix<lower=0>[G, C] mu;   // means of MVN distribution
//     matrix[C, C] L_Sigma[G];    // Cholesky factor of covariance
//     for (g in 1:G) {
//         activities[g] = exposures[g]' * multiplier[g];
//         L_Sigma[g] = diag_pre_multiply(L_sigma[g], L_Omega[g]);
//     }
//     mu = activities * signatures;
// }
// model {
//     for (g in 1:G) {
//         // Priors
//         exposures[g] ~ dirichlet(kappa);
//         multiplier[g] ~ cauchy(0, 1);
//         L_Omega[g] ~ lkj_corr_cholesky(4);
//         L_sigma[g] ~ cauchy(0, 2.5);
// 
//         // Likelihood
//         counts[g] ~ multi_normal_cholesky(mu[g]', L_Sigma[g]);
//     }
// }
