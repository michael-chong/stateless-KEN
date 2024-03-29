data {
  int<lower=1> K; // number of age groups
  int<lower=1> Time; // number of periods

  real<lower=0, upper=1> prop_f;
  
  vector<lower=0>[K] pop_0_prop_data; // initial population 
  
  // age-specific survival
  row_vector<lower=0>[K] surv_mean[Time]; // a priori age-specific survival
  
  // age-specific fertility
  int<lower=0, upper=K+1> n_fert_ages;
  row_vector<lower=0>[n_fert_ages] fert_mean[Time]; 
  int<lower=0, upper=K+1>ages_before_fert;
  int<lower=0, upper=K+1>ages_after_fert;
  
  //vector[K] migr_mean[Time]; // a priori net migration

  int<lower=1> n_obs;
  int<lower=0, upper=Time> t_obs[n_obs]; 
  vector<lower=0>[K] pop_obs[n_obs];

}

transformed data {
  vector[K] log_pop_0_ratios_data;
  
  row_vector[K] logit_surv_mean[Time] = logit(surv_mean);
  row_vector[n_fert_ages] log_fert_mean[Time] = log(fert_mean);
  
  log_pop_0_ratios_data = log(pop_0_prop_data / pop_0_prop_data[10]);
}

parameters {
  
  // total population and proportions
  real log_pop_0_total_true_raw;
  vector[K] log_pop_0_ratios_raw;
  //real<lower=0> sigma_ratios;
  
  // rates
  row_vector[K] logit_surv_raw[Time];
  row_vector[n_fert_ages] log_fert_sub_raw[Time];
  
  real<lower=0> sigma_pop;
  real<lower=0> sigma_fert;
  real<lower=0> sigma_surv;
}

transformed parameters {
  
  real log_pop_0_total_true;
  real<lower=0> pop_0_total_true;
  
  vector[K] pop_0_props;
  
  // transforms of (survival rates)
  row_vector[K] logit_surv[Time];
  row_vector<lower=0, upper=1>[K+1] surv[Time];
  
  // (transforms of) fertility rate
  row_vector[n_fert_ages] log_fert_sub[Time];
  row_vector<lower=0>[n_fert_ages] fert_sub[Time];
  row_vector<lower=0>[K+1] fert[Time];
  
  // leslie matrix and true population
  vector<lower=0>[K] pop_true[Time];
  matrix[1, K] reprod[Time];
  matrix[K, K] leslie_mat[Time];
  matrix[K-1, K] leslie_mat_sub[Time]; 

  // initial population -----------------------------------
  pop_0_props = softmax(log_pop_0_ratios_raw*1+ log_pop_0_ratios_data);
  
  pop_0_total_true = exp(log_pop_0_total_true_raw*.7 + 7);
  pop_true[1] = pop_0_total_true*pop_0_props;

  // leslie matrix construction ----------------------------
  // non-centered parameterization of logit_surv
  for (t in 1:Time) {
    logit_surv[t] = logit_surv_mean[t] + logit_surv_raw[t]*sigma_surv;
  }

  // un-transform survival from logit scale 
  for (t in 1:Time) {
    surv[t] = append_col(inv_logit(logit_surv[t]), 0);
  }
  
  // build fertility rates
  for (t in 1:Time) {
    log_fert_sub[t] = log_fert_mean[t] + log_fert_sub_raw[t]*sigma_fert;
  }
  
  fert_sub = exp(log_fert_sub);

  for (t in 1:Time) {
    fert[t] = append_col(append_col(rep_row_vector(0, ages_before_fert), fert_sub[t]), rep_row_vector(0, ages_after_fert + 1));
  }
  
  for (t in 1:Time) {
  
    // first row of leslie matrix
    reprod[t] = [2.5*surv[t, 1]*prop_f*(fert[t, 1:K] + fert[t, 2:(K+1)] .* surv[t, 2:(K+1)])];
    
    leslie_mat_sub[t] = append_col(diag_matrix(to_vector(surv[t, 2:K])), append_row(rep_vector(0, K-2), surv[t, K+1]));
    leslie_mat[t] = append_row(reprod[t], leslie_mat_sub[t]);

  }  
  
  for (t in 2:Time) {
    pop_true[t] = leslie_mat[t-1] * pop_true[t-1];
  }  
}

model {
  
  log_pop_0_ratios_raw ~ std_normal();
  log_pop_0_total_true_raw ~ std_normal();
  
  if (t_obs[1] > 0) {
    for (i in 1:n_obs) {
      pop_obs[i] ~ normal(pop_true[t_obs[i]], sigma_pop);
    }
  }

  // priors for rates
  for (t in 1:Time) {
    log_fert_sub_raw[t] ~ std_normal();
  }

  for (t in 1:Time) {
    logit_surv_raw[t] ~ std_normal();
  }
  
  sigma_pop ~ normal(0, 10);
  sigma_fert ~ normal(0, 10);
  sigma_surv ~ normal(0, 10);
}

generated quantities {
  real total_pop[Time];
  
  for (t in 1:Time) {
    total_pop[t] = sum(pop_true[t]);
  }
}
