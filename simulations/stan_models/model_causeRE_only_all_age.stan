data {
    int<lower=1> N; // number of observations
    int<lower=1> p; // number of covariates
    matrix[N,p] X; // model matrix for intercepts
    int<lower=1> ncause; // number of causes
    int<lower=1> s; // number of tabulations across which to have correlated cause REs
    int<lower=0> s_index[N]; // indicator of non-cause tabulation index (1:s)
    int<lower=0> c[N]; // indicator of causes (1:ncause)
    vector<lower=0>[N] logpy; // person years in VR
    int<lower=0> y[N]; // data of tabulated VR deaths
    int<lower=0> lkj_hyperpar; // hyperparameter for the LKJ prior on correlations
}
parameters {
    cholesky_factor_corr[ncause] Lcorr; // cholesky factor of the correlations between causes
    vector<lower=0>[ncause] sigma; // standard deviations of the cause REs
    matrix[s,ncause] b_star; // cause-specific RE term (to be transformed for easier sampling)
    vector[p] beta; // log-scale strata-level covariate coefficient parameter
}
transformed parameters {
    vector[N] alpha; // X%*%beta
    vector[N] b; // vector of b's, each cause appended
    matrix[s,ncause] b_mat; // transformed b
    cov_matrix[ncause] D; // covariance of b
    matrix[ncause,ncause] Omega; // correlation matrix for b
    matrix[ncause, ncause] L; // cholesky decomposition of D

    Omega = multiply_lower_tri_self_transpose(Lcorr);
    D = quad_form_diag(Omega, sigma);
    L = cholesky_decompose(D);
    
    alpha = X*beta; // the strata-level mean (fixed effects)
    
    b_mat = b_star * L'; // b in matrix form

    for (i in 1:N) {
        b[i] = b_mat[s_index[i],c[i]];
    }
}
model {
    sigma ~ student_t(3,0,1); // half t prior on the standard deviations (because sigma is lower bounded by 0)
    Lcorr ~ lkj_corr_cholesky(lkj_hyperpar); // LKJ prior recommended by Stan authors
    beta ~ normal(0,5); // if these aren't specified, they have an improper uniform prior
    to_vector(b_star) ~ normal(0,1); // implies b ~ multi_normal(0, D)

    y ~ poisson_log(alpha + b + logpy); // given lambda, the VR deaths are poisson
}
generated quantities {
    vector[N] lambda; // predicted conditional log mortality rate
    vector[N] pred_deaths; // predicted deaths
    vector[N] log_lik; // log likelihood for use with WAIC or PSIS-LOO
    real log_lambda_star; // log marginal all-cause mortality rate
    vector[ncause] log_lambda_tilde; // log marginal cause-specific mortality rate
    // vector[ncause] csmfs; // marginal CSMFs

    lambda = alpha + b; // predicted log mort rate
    pred_deaths = exp(logpy + lambda); // predicted deaths
    for (n in 1:N) {
        log_lik[n] = poisson_log_lpmf(y[n] | alpha[n] + b[n] + logpy[n]);
    }
    log_lambda_tilde = mean(alpha) + (0.5 * sigma .* sigma);
    log_lambda_star = log(sum(exp(log_lambda_tilde)));
    // csmfs = exp(log_lambda_tilde)/exp(log_lambda_star);
}

