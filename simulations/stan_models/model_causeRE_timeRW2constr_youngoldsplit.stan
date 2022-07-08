data {
    int<lower=1> N; // number of observations
    int<lower=1> p; // number of covariates
    matrix[N,p] X; // model matrix for intercepts
    int<lower=1> ncause; // number of causes
    int<lower=1> ncause_young; // number of causes for young age groups
    int<lower=1> ncause_old; // number of causes for old age groups (always a subset of the young causes!!)
    int<lower=0> young_causes[ncause_young]; // indicator of which causes are young for covariance matrix
    int<lower=0> old_causes[ncause_old]; // indicator of which causes are old for covariance matrix
    int<lower=1> s_young; // number of tabulations across which to have correlated cause REs for younger age groups
    int<lower=1> s_old; // number of tabulations across which to have correlated cause REs for older age groups
    int<lower=0> s_index_young[N]; // indicator of non-cause tabulation index (1:s_young)
    int<lower=0> s_index_old[N]; // indicator of non-cause tabulation index (1:s_old)
    int<lower=0> c_young[N]; // indicator of young causes (1:ncause_young)
    int<lower=0> c_old[N]; // indicator of old causes (1:ncause_old)
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
    matrix[s_young,ncause_young] b_star_young; // cause-specific RE term for young (to be transformed for easier sampling)
    matrix[s_old,ncause_old] b_star_old; // cause-specific RE term for old (to be transformed for easier sampling)
    vector[p] beta; // log-scale strata-level covariate coefficient parameter
    matrix[n_rw_index,nyear] gamma_mat; // log-scale time random walk parameter
    real<lower=0> sigma_gamma; // time random walk std dev
}
transformed parameters {
    vector[N] alpha; // X%*%beta
    vector[N] b; // vector of b's, each cause appended
    vector[N] gamma_vec; // long vector of time random walk parameters
    matrix[s_young,ncause_young] b_mat_young; // transformed b for young causes
    matrix[s_old,ncause_old] b_mat_old; // transformed b for old causes
    cov_matrix[ncause] D; // covariance of b
    matrix[ncause,ncause] Omega; // correlation matrix for b
    matrix[ncause, ncause] L; // cholesky decomposition of D
    
    Omega = multiply_lower_tri_self_transpose(Lcorr);
    D = quad_form_diag(Omega, sigma);
    L = cholesky_decompose(D);
    
    alpha = X*beta; // the strata-level mean (fixed effects)
    
    b_mat_young = b_star_young * L[young_causes,young_causes]'; // b in matrix form for young causes 
    b_mat_old = b_star_old * L[old_causes,old_causes]'; // b in matrix form for old causes 

    for (i in 1:N) {
        gamma_vec[i] = gamma_mat[rw_index[i],k[i]];
        if (s_index_young[i] != 0) {
            b[i] = b_mat_young[s_index_young[i],c_young[i]];
        } else {
            b[i] = b_mat_old[s_index_old[i],c_old[i]];
        }
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
    Lcorr ~ lkj_corr_cholesky(lkj_hyperpar); // LKJ prior on correlation matrix recommended by Stan  authors
    sigma_gamma ~ student_t(3,0,1); // similar half t as above
    beta ~ normal(0,10); // if these aren't specified then it means we have an improper uniform prior on all betas

    to_vector(b_star_young) ~ normal(0,1); // implies b ~ multi_normal(0, D)
    to_vector(b_star_old) ~ normal(0,1); // implies b ~ multi_normal(0, D)
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

