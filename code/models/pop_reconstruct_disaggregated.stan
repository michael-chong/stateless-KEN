
data {
  int<lower=1> K; // number of age groups
  int<lower=1> Time; // number of periods
  
  real<lower=0, upper=1> prop_f; // proportion female at birth
  
  // age-specific survival
  array[Time] row_vector<lower=0>[K] surv_mean; // a priori age-specific survival
  
  // age-specific fertility
  int<lower=0, upper=K + 1> n_fert_ages;
  array[Time] row_vector<lower=0>[n_fert_ages] fert_mean;
  int<lower=0, upper=K + 1> ages_before_fert;
  int<lower=0, upper=K + 1> ages_after_fert;
  
  // population age proportions from origin country
  array[Time] vector<lower=0>[K] init_mig_pop_props_data;
  
  int<lower=1> n_obs;
  array[n_obs] int<lower=0, upper=Time> t_obs;
  array[n_obs] vector<lower=0>[K] pop_native_obs;
  array[n_obs, Time] vector<lower=0>[K] pop_mig_obs;
  
  // hyperprior
  real<lower=0> fert_hyper;
  real<lower=0> surv_hyper;
  real<lower=0> sd_init_ratios;
}
transformed data {
  array[Time] vector[K] log_init_mig_pop_ratios_data;
  
  array[Time] row_vector[K] logit_surv_mean = logit(surv_mean);
  array[Time] row_vector[n_fert_ages] log_fert_mean = log(fert_mean);
  
  for (t in 1 : Time) {
    log_init_mig_pop_ratios_data[t] = log(init_mig_pop_props_data[t]
                                          / init_mig_pop_props_data[t, 10]);
  }

}
parameters {
  // total population and proportions
  vector[K] log_pop_mig_init_ratios;
  array[Time] real log_pop_mig_init_total;
  
  //real<lower=0> sigma_ratios;
  
  // rates
  array[Time] row_vector[K] logit_surv_raw;
  array[Time] row_vector[n_fert_ages] log_fert_sub_raw;
  
  //real<lower=0> sigma_pop;
  real<lower=0> sigma_fert;
  real<lower=0> sigma_surv;
}
transformed parameters {
  // transforms of (survival rates)
  array[Time] row_vector[K] logit_surv;
  array[Time] row_vector<lower=0, upper=1>[K + 1] surv;
  
  // (transforms of) fertility rate
  array[Time] row_vector[n_fert_ages] log_fert_sub;
  array[Time] row_vector<lower=0>[n_fert_ages] fert_sub;
  array[Time] row_vector<lower=0>[K + 1] fert;
  
  // initial population
  array[Time] real<lower=0> pop_mig_init_total;
  array[Time] vector<lower=0>[K] pop_mig_init_props;
  
  // leslie matrix and true population
  array[Time] vector<lower=0>[K] pop_native_true;
  array[Time, Time] vector[K] pop_mig_true; // [mig year, actual year, age]
  array[Time, Time] real<lower=0> mig_births;
  array[Time] row_vector[K] reprod;
  array[Time] matrix[K, K] leslie_mat;
  array[Time] matrix[K, K] leslie_mat_nb;
  array[Time] matrix[K - 1, K] leslie_mat_sub;
  
  // leslie matrix construction ----------------------------
  // non-centered parameterization of logit_surv
  for (t in 1 : Time) {
    logit_surv[t] = logit_surv_mean[t] + logit_surv_raw[t] * sigma_surv;
  }
  
  // un-transform survival from logit scale 
  for (t in 1 : Time) {
    surv[t] = append_col(inv_logit(logit_surv[t]), 0);
  }
  
  // build fertility rates
  for (t in 1 : Time) {
    log_fert_sub[t] = log_fert_mean[t] + log_fert_sub_raw[t] * sigma_fert;
  }
  
  fert_sub = exp(log_fert_sub);
  
  for (t in 1 : Time) {
    fert[t] = append_col(append_col(rep_row_vector(0, ages_before_fert),
                                    fert_sub[t]),
                         rep_row_vector(0, ages_after_fert + 1));
  }
  
  for (t in 1 : Time) {
    // first row of leslie matrix
    reprod[t] = to_row_vector([2.5 * surv[t, 1] * prop_f
                               * (fert[t, 1 : K]
                                  + fert[t, 2 : (K + 1)]
                                    .* surv[t, 2 : (K + 1)])]);
    
    // all other rows of leslie matrix
    leslie_mat_sub[t] = append_col(diag_matrix(to_vector(surv[t, 2 : K])),
                                   append_row(rep_vector(0, K - 2),
                                              surv[t, K + 1]));
    
    // entire leslie matrix
    leslie_mat[t] = append_row(reprod[t], leslie_mat_sub[t]);
    
    // leslie matrix without births
    leslie_mat_nb[t] = append_row(rep_row_vector(0, K), leslie_mat_sub[t]);
  }
  
  // initialize migration populations
  pop_mig_init_total = exp(log_pop_mig_init_total);
  // set the most recent migration wave to 0 
  for (t in (t_obs[n_obs]+1):Time) {
    pop_mig_init_total[t] = 0;
  }

  pop_mig_true = rep_array(rep_vector(0, K), Time, Time);
  
  for (mig_year in 1 : Time) {
    pop_mig_init_props[mig_year] = softmax(log_pop_mig_init_ratios[mig_year]
                                           * sd_init_ratios
                                           + log_init_mig_pop_ratios_data[mig_year]);
    
    pop_mig_true[mig_year, mig_year] = pop_mig_init_total[mig_year]
                                       * pop_mig_init_props[mig_year];
  }
  
  mig_births = rep_array(0, Time, Time);
  
  // aging and deaths for migrant population
  for (mig_year in 1 : (Time - 1)) {
    for (t in (mig_year + 1) : Time) {
      pop_mig_true[mig_year, t] = leslie_mat_nb[t-1]
                                  * pop_mig_true[mig_year, t - 1];
      mig_births[t, mig_year] = dot_product(reprod[t-1],
                                            pop_mig_true[mig_year, t - 1]);
    }
  }
  
  pop_native_true[1] = rep_vector(0, K);
  
  for (t in 2 : Time) {
    pop_native_true[t] = leslie_mat[t-1] * pop_native_true[t - 1]
                         + append_row(sum(mig_births[t]),
                                      rep_vector(0, K - 1));
  }
}
model {
  if (t_obs[1] > 0) {
    for (i in 1 : n_obs) {
      // migrant populations
      for (mig_year in 1 : Time) {
        pop_mig_obs[i, mig_year] ~ normal(pop_mig_true[mig_year, t_obs[i]],
                                          1);
      }
      
      // native-born population
      pop_native_obs[i] ~ normal(pop_native_true[t_obs[i]], 1);
    }
  }
  
  // priors for rates
  for (t in 1 : Time) {
    log_fert_sub_raw[t] ~ std_normal();
  }
  
  for (t in 1 : Time) {
    logit_surv_raw[t] ~ std_normal();
  }
  
  // initial population
  log_pop_mig_init_total ~ normal(4, 1);
  log_pop_mig_init_ratios ~ std_normal();
  
  //sigma_pop ~ normal(0, 100);
  sigma_fert ~ normal(0, fert_hyper);
  sigma_surv ~ normal(0, surv_hyper);
}
generated quantities {
  array[Time] real total_mig_pop;
  array[Time] real total_native_pop;
  array[Time] real total_pop;
  
  for (t in 1 : Time) {

    total_native_pop[t] = sum(pop_native_true[t]);
    
    total_mig_pop[t] = 0;
    for (mig_year in 1 : Time) {
      total_mig_pop[t] = total_mig_pop[t] + sum(pop_mig_true[mig_year, t]);
    }
    
    total_pop[t] = total_native_pop[t] + total_mig_pop[t];
    
  }
  
}


