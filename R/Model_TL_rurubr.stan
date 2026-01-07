data {
  int<lower=1> N;             // n individuals
  int<lower=1> H;             // n households
  int<lower=1> V;             // n villages
  array[N] int<lower=1, upper=H> hh; // household index for each individual (1..H)
  array[H] int<lower=1, upper=V> vv; // village index for each household (1..V)
  vector<lower=0>[N] age;     // vector with age for each individual
  array[N] int<lower=0, upper=1> y;  // serological result (1 = pos, 0 = neg)
  real national_mu;
  real national_sd;
  real mean_sigma_n;
  real mean_sigma_v;
  real sd_sigma_n;
  real sd_sigma_v;
  array[V] real X_v; //1 if urban, 0 if rural
}

parameters {
  real<lower=0> sigma_n;       // between-village SD (on logit scale)
  real<lower=0> sigma_v;       // between-household SD (on logit scale)
  vector[V] delta_v;           // village-level random effect
  vector[H] delta_h;           // household-level random effect
  real<lower=0,upper=1> se;      // sensitivity
  real<lower=0,upper=1> sp;      // specificity
  real beta;  //urban effect
}

transformed parameters {
  vector[V] mu_v;              // logit(FOI) per village
  vector[V] lambda_v;   // FOI per village (probability)
  vector[H] logit_lambda_h;    // logit(FOI) per household
  vector[H] lambda_h;   // FOI per household
  vector[N] lambda_i;          // FOI per individual
  vector[N] z;                 // true seroprevalence
  vector[N] p_obs;             // observed seroprevalence (with se/sp)

  // --- VILLAGE LEVEL ---
  for (v in 1:V) {
    mu_v[v] = national_mu + delta_v[v] + beta * X_v[v];           // logit of FOI
    lambda_v[v] = inv_logit(mu_v[v]);         
  }

  // --- HOUSEHOLD LEVEL ---
  for (h in 1:H) {
    logit_lambda_h[h] = mu_v[vv[h]] + delta_h[h]; // logit of household FOI
    lambda_h[h] = inv_logit(logit_lambda_h[h]);  
  }

  // --- INDIVIDUAL LEVEL ---
  for (i in 1:N) {
    lambda_i[i] = lambda_h[hh[i]];
    z[i] = 1 - exp(- lambda_i[i] * age[i]);
    p_obs[i] = se * z[i] + (1 - sp) * (1 - z[i]);
  }
}

model {
  // --- PRIORS ---
  sigma_n ~ normal(mean_sigma_n, sd_sigma_n);
  sigma_v ~ normal(mean_sigma_v, sd_sigma_v);
  beta ~ normal(0, 1);
  delta_v ~ normal(0, sigma_n);
  delta_h ~ normal(0, sigma_v);
  se ~ normal(0.77, 0.1); //sensitivity around different sd 0.05
  sp ~ normal(0.93, 0.1); //sensitivity around different sd 0.05
  
  // --- LIKELIHOOD ---
  for (i in 1:N)
   target += bernoulli_lpmf(y[i] | p_obs[i]);
}

generated quantities {
  vector[N] log_lik;
  for (i in 1:N)
    log_lik[i] = bernoulli_lpmf(y[i] | p_obs[i]);
}

