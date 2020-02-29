data {
    int<lower=1> N; // number of observations
    int<lower=1> p; // number of covariates
    matrix[N,p] X; // model matrix for intercepts
    int<lower=0> n_rw_index; // number of tabulations for the RW on time to go over
    int<lower=1> nyear; // number of years
    int<lower=0> k[N]; // indicator of years (1:nyear)
    int<lower=0> rw_index[N]; // indicator of tabulations across which time RW happens (1:n_rw_index)
    vector<lower=0>[N] logpy; // person years in VR
    int<lower=0> y[N]; // data of tabulated VR deaths
}
parameters {
    vector[p] beta; // log-scale strata-level covariate coefficient parameter
    matrix[n_rw_index,nyear] gamma_mat; // log-scale time random walk parameter
    real<lower=0> sigma_gamma; // time random walk std dev
}
transformed parameters {
    vector[N] alpha; // X%*%beta
    vector[N] gamma_vec; // long vector of time random walk parameters
    
    alpha = X*beta; // the strata-level mean (fixed effects)
    
    for (i in 1:N) {
        gamma_vec[i] = gamma_mat[rw_index[i],k[i]];
    }
}
model {
    y ~ poisson_log(alpha + gamma_vec + logpy); // given lambda, the VR deaths are poisson
    for (l in 1:n_rw_index) {
        for (m in 3:(nyear)) {
            gamma_mat[l,m] ~ normal(2*gamma_mat[l,m-1] - gamma_mat[l,m-2],sigma_gamma); // RW2 conditional; the first two have improper uniform dist'ns
        }
        sum(gamma_mat[l,]) ~ normal(0, 0.00001 * nyear);  // equivalent to mean(gamma_mat[l,]) ~ normal(0,0.00001)
    }
    
    sigma_gamma ~ student_t(3,0,1); // half t distribution on std dev
    beta ~ normal(0,5); // if these aren't specified then it means we have an improper uniform prior on all betas
}
generated quantities {
    vector[N] lambda; // predicted conditional log mortality rate
    vector[N] mu; // predicted conditional mortality rate
    vector[N] pred_deaths; // predicted deaths
    vector[N] log_lik; // log likelihood for use with WAIC or PSIS-LOO

    lambda = alpha + gamma_vec; // predicted log mort rate
    mu = exp(lambda);
    pred_deaths = mu .* exp(logpy); // predicted deaths
    for (n in 1:N) {
        log_lik[n] = poisson_log_lpmf(y[n] | alpha[n] + gamma_vec[n] + logpy[n]);
    }
}

