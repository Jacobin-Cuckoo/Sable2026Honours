functions { }

data {
  int<lower=1> N;                 // number of subjects
  int<lower=1> n_obs;             // longitudinal observations
  int<lower=1> P;                 // fixed effects
  int<lower=1> Q;                 // random effects
  
  matrix[n_obs,P] X;
  
  matrix[n_obs,Q] Z;
  
  vector[n_obs] y;
  
  array[n_obs] int id;
  
  vector[N] surv_time;
  
  array[N] int event;
  
  matrix[N,P] W;
  
}

parameters {
  
  vector[P] beta;
  
  vector[P] gamma;
  
  vector[Q] b[N];
  
  cholesky_factor_corr[Q] L;
  
  vector<lower=0>[Q] tau;
  
  real<lower=0> sigma;
  
  real alpha;
  
}

transformed parameters {
  
  matrix[Q,Q] D;
  
  D =
    diag_pre_multiply(tau,L)
  *
    diag_pre_multiply(tau,L)';

}

model {

  //---------------------------------------------------
  // Priors
  //---------------------------------------------------

  beta ~ normal(0,5);

  gamma ~ normal(0,5);

  alpha ~ normal(0,2);

  tau ~ cauchy(0,2);

  sigma ~ cauchy(0,2);

  L ~ lkj_corr_cholesky(2);

  //---------------------------------------------------
  // Random effects
  //---------------------------------------------------

  for(i in 1:N)
      b[i] ~ multi_normal_cholesky(rep_vector(0,Q),
                                   diag_pre_multiply(tau,L));

  //---------------------------------------------------
  // Longitudinal model
  //---------------------------------------------------

  for(n in 1:n_obs){

      real mu;

      mu =
      X[n]*beta
      +
      Z[n]*b[id[n]];

      y[n]
      ~
      normal(mu,sigma);

  }

  //---------------------------------------------------
  // Cox likelihood
  //---------------------------------------------------

  for(i in 1:N){

      real eta;

      eta =
      W[i]*gamma
      +
      alpha*
      (
      W[i]*beta
      +
      W[i]*b[i]
      );

      if(event[i]==1){

         target += eta;

      }

  }

}