// https://discourse.mc-stan.org/t/multivariate-likelihoods/1080/7




// // https://groups.google.com/d/msg/stan-users/nk03n8hX4Us/-pTN5IGjCQAJ
// functions {
//     //#include "common_functions.stan"
//     matrix array_to_matrix(vector[] x) {
//         // Assume x doesn't have 0 rows or columns
//         matrix[size(x), rows(x[1])] y;
//         for (m in 1:size(x))
//             y[m] = x[m]';
//         return y;
//     }
// }
// data {
//     int<lower=1> C;                // number of mutation categories
//     int<lower=1> S;                // number of mutational signatures
//     int<lower=1> G;                // number of genomes
//     matrix[S, C] signatures;       // matrix of signatures (rows) to be fitted
//     matrix<lower=0>[G, C] counts;  // matrix of counts per category (columns) per genome sample (rows)
//     vector<lower=0>[S] kappa;      // prior on exposures (i.e. mixing proportions of signatures)
// }
// parameters {
//     simplex[S] exposures[G];             // parameters of interest (mixing proportions)
//     real<lower=1> nu[G];                 // degrees of freedom
//     cholesky_factor_corr[C] L_Omega[G];  // Cholesky factor of correlation (Stan manual, S9.15)
//     matrix<lower=0>[G, C] L_sigma;       // standard deviations
//     matrix[G, C] z;                      // for stochastic representation of MVT
//     real<lower=0> u[G];                  // idem
// }
// transformed parameters {
//     matrix[G, C] x;                        // for stochastic representation of MVT; x[g] are distributed MVT (multi_student_t)
//     matrix<lower=0, upper=1>[G, C] mu;     // idem
//     matrix[C, C] L_Sigma[G];               // Cholesky factor of covariance
//     matrix<lower=0, upper=1>[G, C] probs;  // mutation type probabilities
//     probs = array_to_matrix(exposures) * signatures;
//     // (array_to_matrix is defined in common_functions.stan and is not in base Stan)
//     for (g in 1:G) {
//         L_Sigma[g] = diag_pre_multiply(L_sigma[g], L_Omega[g]);
//         mu[g] = probs[g] * sum(counts[g]);
//         x[g] = to_row_vector(to_vector(mu[g]) + sqrt(nu[g] / u[g]) * (L_Sigma[g] * to_vector(z[g])));  // distributed MVT
//     }
// }
// model {
//     nu ~ gamma(2, 0.1);
//     for (g in 1:G) {
//         exposures[g] ~ dirichlet(kappa);
//         target += normal_lpdf(z[g] | 0, 1);
//         target += chi_square_lpdf(u[g] | nu[g]);
//     }
// }
