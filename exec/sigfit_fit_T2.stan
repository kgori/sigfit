functions {
//    #include "common_functions.stan"
}
data {
    int<lower=1> C;               // number of mutation categories
    int<lower=1> S;               // number of mutational signatures
    int<lower=1> G;               // number of genomes
    matrix[G, C] counts;          // observed mutation counts (row per genome)
    matrix[S, C] signatures;      // signatures to fit (row per signature)
    vector<lower=0>[S] kappa;     // prior on exposures (mixing proportions)
    matrix[G, C] opps;            // mutational opportunities (row per genome)
}
parameters {
    simplex[S] exposures[G];      // signature exposures (row per genome)
    // Auxiliary parameters for reparameterising half-Cauchy parameters
    vector<lower=0, upper=pi()/2>[G] multiplier_unif;
    vector<lower=0, upper=pi()/2>[G] sigma_unif;
    vector<lower=0, upper=pi()/2>[C] nu_unif;
}
transformed parameters {
    matrix<lower=0>[G, S] activities;  // signature activities (row per genome)
    matrix[G, C] expected_counts;      // normal means (row per genome)
    vector<lower=0>[G] multiplier;     // exposure multipliers
    vector<lower=0>[G] sigma;          // standard deviations
    vector<lower=0>[C] nu;             // degrees of freedom
    // Reparameterisation of multiplier, sigma [Half-Cauchy(0, 2.5)]
    multiplier = 2.5 * tan(multiplier_unif);
    sigma = 2.5 * tan(sigma_unif);
    nu = 2.5 * tan(nu_unif);
    for (g in 1:G) {
        activities[g] = exposures[g]' * multiplier[g];
    }
    expected_counts = activities * signatures .* opps;
}
model {
    //multiplier ~ cauchy(0, 2.5);  // multiplier priors
    //sigma ~ cauchy(0, 2.5);       // std. dev. priors
    for (g in 1:G) {
        exposures[g] ~ dirichlet(kappa);                             // exposure priors
        for (c in 1:C) {
            counts[g, c] ~ student_t(nu[c], expected_counts[g, c], sigma[g]);  // t likelihood
        }
    }
}
