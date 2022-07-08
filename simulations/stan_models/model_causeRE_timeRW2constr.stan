data {
    int<lower=1> N; // number of observations (ncause*nage*nyear*nstrata)
    int<lower=1> p; // number of covariates
    matrix[N,p] X; // model matrix for fixed effects
    int<lower=1> ncause; // number of causes
    int<lower=1> s; // number of tabulations across which to have correlated cause REs
    int<lower=0> s_index[N]; // indicator of non-cause tabulation index (1:s)
    int<lower=0> c[N]; // indicator of causes (1:ncause)
    int<lower=0> n_rw_index; // number of tabulations for the RW on time to go over
    int<lower=1> nyear; // number of years
    int<lower=0> k[N]; // indicator of years (1:nyear)
    int<lower=0> rw_index[N]; // indicator of tabulations across which time RW happens (1:n_rw_index)
    vector<lower=0>[N] logpy; // person years in VR
    int<lower=0> y[N]; // data of tabulated VR deaths
    int<lower=0> lkj_hyperpar; // hyperparameter for the LKJ prior
}
parameters {
    cholesky_factor_corr[ncause] Lcorr; // cholesky factor of the correlations between causes
    vector<lower=0>[ncause] sigma; // standard deviations of the cause REs
    matrix[s,ncause] b_star; // cause-specific RE term (to be transformed for easier sampling)
    vector[p] beta; // log-scale strata-level covariate coefficient parameter
    matrix[n_rw_index,nyear] gamma_mat; // log-scale time random walk parameter
    real<lower=0> sigma_gamma; // time random walk std dev
}
transformed parameters {
    vector[N] alpha; // X%*%beta
    vector[N] b; // vector of b's, each cause appended
    vector[N] gamma_vec; // long vector of time random walk parameters
    matrix[s,ncause] b_mat; // transformed b as a matrix
    cov_matrix[ncause] D; // covariance of b (reffered to as Sigma, generally)
    matrix[ncause,ncause] Omega; // correlation matrix for b
    matrix[ncause, ncause] L; // cholesky decomposition of D
    
    Omega = multiply_lower_tri_self_transpose(Lcorr);
    D = quad_form_diag(Omega, sigma);
    L = cholesky_decompose(D);
    
    alpha = X*beta; // the strata-level mean (fixed effects)
    
    b_mat = b_star * L'; // b in matrix form

    for (i in 1:N) {
        gamma_vec[i] = gamma_mat[rw_index[i],k[i]];
        b[i] = b_mat[s_index[i],c[i]];
    }
}
model {
    y ~ poisson_log(alpha + gamma_vec + b + logpy); // given lambda, the VR deaths are poisson
    for (l in 1:n_rw_index) {
        for (m in 3:(nyear)) {
            gamma_mat[l,m] ~ normal(2*gamma_mat[l,m-1] - gamma_mat[l,m-2],sigma_gamma); // RW2 conditional; the first two have improper uniform dist'ns
        }
        sum(gamma_mat[l,]) ~ normal(0, 0.00001 * nyear);  // equivalent to mean(gamma_mat[l,]) ~ normal(0,0.00001)
    }
    
    sigma ~ student_t(3,0,1); // leads to a half t prior on the standard deviations (because sigma is lower bounded by 0)
    Lcorr ~ lkj_corr_cholesky(lkj_hyperpar); // LKJ prior recommended by Stan authors
    sigma_gamma ~ student_t(3,0,1); // half t, similar to above
    beta ~ normal(0,5); // if these aren't specified then it means we have an improper uniform prior on all betas

    to_vector(b_star) ~ normal(0,1); // implies b ~ multi_normal(0, D)
}
generated quantities {
    vector[N] lambda; // predicted conditional log mortality rate
    vector[N] mu; // predicted conditional mortality rate
    vector[N] pred_deaths; // predicted deaths
    vector[N] log_lik; // log likelihood for use with WAIC or PSIS-LOO

    lambda = alpha + b + gamma_vec; // predicted log mort rate
    mu = exp(lambda);
    pred_deaths = mu .* exp(logpy); // predicted deaths
    for (n in 1:N) {
        log_lik[n] = poisson_log_lpmf(y[n] | alpha[n] + b[n] + gamma_vec[n] + logpy[n]);
    }
}

