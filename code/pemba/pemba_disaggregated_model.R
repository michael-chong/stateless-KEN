library(tidyverse)
library(here)
library(readxl)
library(cmdstanr)
library(readxl)
library(posterior)

pemba_raw <- read_excel(here("data/unhcr/kenya_pemba/unhcr_ken_2016_pemba-anonymized.xls"))

pemba_tribes <- read_csv(here("data/unhcr/kenya_pemba/kenya_tribes.csv"))

pemba <- pemba_raw %>%
  filter(!tribe %in% pemba_tribes$`Tribe name`) %>%
  filter(!is.na(age_member)) 

all_pop <- read_excel(here("data/wpp/WPP2019_POP_F07_3_POPULATION_BY_AGE_FEMALE.xlsx"), skip = 16) %>%
  rename(
    region = `Region, subregion, country or area *`,
    period = `Reference date (as of 1 July)`
  ) %>%
  select(-Index, -Variant, -Notes, -`Country code`, -`Type`, -`Parent code`) %>% 
  pivot_longer(names_to = "age", values_to = "pop", cols = -c("region", "period")) %>%
  mutate(age = as.numeric(str_extract(age, ".+?(?=[-+])")), pop = as.numeric(pop)) 

# prepare fertility and mortality rates from WPP
all_fert <- read_excel(here("data/wpp/WPP2019_FERT_F07_AGE_SPECIFIC_FERTILITY.xlsx"), skip = 16) %>%
  rename(region = `Region, subregion, country or area *`, period = Period) %>% 
  select(-Index, -Variant, -Notes, -`Country code`, -`Parent code`, -`Type`) %>% 
  mutate(year = as.numeric(str_extract(period, ".+?(?=-)"))) %>% # match everything before hyphen or plus
  pivot_longer(names_to = "age", values_to = "fert", cols = -c("region", "period", "year")) %>% 
  mutate(age = as.numeric(str_extract(age, ".+?(?=-)")), fert = as.numeric(fert)/1000) %>%
  select(-period)

all_surv <- read_excel(here("data/wpp/WPP2019_MORT_F17_3_ABRIDGED_LIFE_TABLE_FEMALE.xlsx"), skip = 16) %>%
  transmute(
    region = `Region, subregion, country or area *`,
    age = `Age (x)`,
    period = Period,
    Lx = as.numeric(`Number of person-years lived L(x,n)`)
  ) %>%
  group_by(region, period) %>%
  mutate(surv = Lx / lag(Lx)) %>%
  mutate(year = as.numeric(str_extract(period, ".+?(?=-)")))  %>%
  # adjustments to age groups
  mutate(age = ifelse(age == 1, 0, age)) %>%
  group_by(age) %>%
  group_by(region, year, age) %>%
  summarise(Lx = sum(Lx)) %>%
  mutate(surv = lead(Lx) / Lx) %>%
  ungroup()

kenya_surv <- all_surv %>%
  filter(region == "Kenya") %>%
  filter(!is.na(surv))

kenya_fert <- all_fert %>%
  filter(region == "Kenya") %>%
  right_join(select(kenya_surv, -surv, -Lx)) %>%
  arrange(region, year, age) %>%
  mutate(fert = ifelse(is.na(fert), 0, fert))

n_ages <- length(unique(kenya_fert$age))
n_periods <- length(unique(kenya_fert$year))
n_obs <- 1 # number of population observations (min. 1 for model to work)
# note: set t_obs = as.array(0) below if no actual observations

years <- unique(kenya_fert$year)
years <- c(years, max(years + 5))
ages <- c(unique(kenya_surv$age), Inf)

pemba_mig <- pemba %>%
  filter(!is.na(mig) | mig != "Non applicable") %>%
  mutate(mig_year = as.numeric(mig)) %>%
  filter(mig_year > min(years)) %>%
  count(age_member, mig_year) %>%
  mutate(age_member = str_replace(age_member, "Below 15", "0-15")) %>%
  separate(age_member, into = c("age_start_range", "age_end_range")) %>%
  mutate(
    age_start_range = as.numeric(age_start_range),
    age_end_range = as.numeric(age_end_range)
  ) %>%
  mutate(range_length = age_end_range - age_start_range) %>%
  group_by(age_start_range, mig_year) %>%
  summarize(
    age = age_start_range:(age_end_range-1),
    n = n/range_length
  ) %>%
  arrange(age) %>%
  mutate(
    mig_year = cut(mig_year, breaks = years, labels = head(years, -1), include.lowest = T, right = F),
    age = cut(age, breaks = ages, labels = head(ages, -1), include.lowest = T, right = F)
  ) %>%
  group_by(mig_year, age) %>%
  summarize(n = sum(n)) %>%
  mutate(mig_year = as.numeric(as.character(mig_year)), age = as.numeric(as.character(age))) %>%
  ungroup() %>%
  complete(mig_year = head(years, -1), age = head(ages, -1), fill = list(n = 0)) 

pop_mig_obs <- matrix(
  data = pemba_mig$n,
  nrow = n_periods,
  ncol = n_ages,
  byrow = TRUE
) %>% array(dim = c(1, dim = c(n_periods, n_ages)))

pemba_native <- pemba %>%
  filter(
    cbirth %in% c("Kenya", "Pemba"),
    !is.na(mig) | mig == "Non applicable"
  ) %>%
  count(age_member) %>%
  mutate(age_member = str_replace(age_member, "Below 15", "0-15")) %>%
  separate(age_member, into = c("age_start_range", "age_end_range")) %>%
  mutate(
    age_start_range = as.numeric(age_start_range),
    age_end_range = as.numeric(age_end_range)
  ) %>%
  mutate(range_length = age_end_range - age_start_range) %>%
  group_by(age_start_range) %>%
  summarize(
    age = age_start_range:(age_end_range-1),
    n = n/range_length
  ) %>%
  ungroup() %>%
  mutate(
    age = cut(age, breaks = ages, labels = head(ages, -1), include.lowest = T, right = F)
  ) %>%
  mutate(age = as.numeric(as.character(age))) %>%
  group_by(age) %>%
  summarize(n = sum(n)) %>%
  complete(age = head(ages, -1), fill = list(n = 0))

pop_native_obs <- matrix(
  data = pemba_native$n,
  nrow = 1,
  ncol = n_ages
)

TZA_props <- all_pop %>%
  filter(region == "United Republic of Tanzania") %>%
  group_by(period) %>%
  mutate(pop_prop = pop/sum(pop)) %>%
  ungroup() %>%
  filter(period %in% head(years, -1), age %in% ages) %>%
  arrange(period, age)

mig_init_props <- matrix(
  TZA_props$pop_prop,
  nrow = n_periods,
  ncol = n_ages,
  byrow = TRUE
)

# put rates into matrix format for stan
kenya_surv_mat <- matrix(kenya_surv$surv, nrow = n_periods, ncol = n_ages, byrow = TRUE) 
kenya_fert_mat <- matrix(kenya_fert$fert, nrow = n_periods, ncol = n_ages, byrow = TRUE) 

# some additional objects to help construct fertility matrix in Stan
fert_ages_bool <- apply(kenya_fert_mat, 2, function(x) any(x > 0))
ages_before_fert = which(fert_ages_bool)[1] - 1
ages_after_fert = n_ages - last(which(fert_ages_bool))

# create mappings between indices and actual values
age_mapping <- tibble(
  number = 1:(n_ages + 1),
  value = c(unique(kenya_fert$age), max(kenya_fert$age) + 5)
)

year_mapping <- data.frame(
  number = 1:n_periods,
  value = unique(kenya_fert$year)
)


# put data in list to pass to Stan
data_list <- list(
  K = n_ages,
  Time = n_periods,
  prop_f = 1/(1 + 1.05),
  # survival rates
  surv_mean = kenya_surv_mat,
  # fertility rates
  fert_mean = kenya_fert_mat[, fert_ages_bool],
  n_fert_ages = sum(fert_ages_bool),
  ages_before_fert = ages_before_fert,
  ages_after_fert = ages_after_fert,
  # population observation
  t_obs = as.array(14),
  n_obs = n_obs,
  pop_mig_obs = pop_mig_obs,
  pop_native_obs = pop_native_obs,
  # for migrant population prior
  init_mig_pop_props_data = mig_init_props,
  # for fert/surv hyper priors
  fert_hyper = .05, # .01
  surv_hyper = .07, # .01 # .07
  sd_init_ratios = 5
)

mod <- cmdstan_model(here("code/models/pop_reconstruct_disaggregated.stan"))

fit <- mod$sample(
  data = data_list,
  seed = 111,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  refresh = 100,
  max_treedepth = 13
  #adapt_delta = 0.9
)

# calculate summary statistics
q_probs <- c(0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975)

summary_table <- fit$summary(
  variables = NULL,
  mean, sd, rhat, ~quantile2(.x, probs = q_probs, na.rm = TRUE), ess_bulk, ess_tail
)

# save results
to_save <- list(
  data = data_list,
  table = summary_table,
  mappings = list(year = year_mapping, age = age_mapping)
)

write_rds(to_save, here("output/intermediate/pemba_disaggregated_model.rds"))
