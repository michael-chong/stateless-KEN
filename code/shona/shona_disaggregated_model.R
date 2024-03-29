library(tidyverse)
library(here)
library(cmdstanr)
library(readxl)
library(posterior)

# Prepare Shona data
load(here("data/shona2019_toronto.rdata"))

wpp_pop_file <- here("data/wpp/WPP2019_POP_F07_3_POPULATION_BY_AGE_FEMALE.xlsx")
all_pop <- bind_rows(
    read_excel(wpp_pop_file, skip = 16, sheet = 1),
    read_excel(wpp_pop_file, skip = 16, sheet = 2)
  ) %>%
  rename(
    region = `Region, subregion, country or area *`,
    period = `Reference date (as of 1 July)`
  ) %>%
  select(-Index, -Variant, -Notes, -`Country code`, -`Type`, -`Parent code`) %>% 
  pivot_longer(names_to = "age", values_to = "pop", cols = -c("region", "period")) %>%
  mutate(age = as.numeric(str_extract(age, ".+?(?=[-+])")), pop = as.numeric(pop)) %>%
  filter(period <= 2020)

# prepare fertility and mortality rates from WPP
wpp_fert_file <- here("data/wpp/WPP2019_FERT_F07_AGE_SPECIFIC_FERTILITY.xlsx")
all_fert <- bind_rows(
    read_excel(wpp_fert_file, skip = 16, sheet = 1),
    read_excel(wpp_fert_file, skip = 16, sheet = 2)
  ) %>%
  rename(region = `Region, subregion, country or area *`, period = Period) %>% 
  select(-Index, -Variant, -Notes, -`Country code`, -`Parent code`, -`Type`) %>% 
  mutate(year = as.numeric(str_extract(period, ".+?(?=-)"))) %>% # match everything before hyphen or plus
  pivot_longer(names_to = "age", values_to = "fert", cols = -c("region", "period", "year")) %>% 
  mutate(age = as.numeric(str_extract(age, ".+?(?=-)")), fert = as.numeric(fert)/1000) %>%
  select(-period) |>
  filter(year <= 2025)

wpp_life_table_file <- here("data/wpp/WPP2019_MORT_F17_3_ABRIDGED_LIFE_TABLE_FEMALE.xlsx")
all_surv <- bind_rows(
    read_excel(wpp_life_table_file, skip = 16, sheet = 1),
    read_excel(wpp_life_table_file, skip = 16, sheet = 2)
  ) %>%
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
  ungroup() |>
  filter(year <= 2025)


kenya_surv <- all_surv %>%
  filter(region == "Kenya") %>%
  filter(!is.na(surv))

kenya_fert <- all_fert %>%
  filter(region == "Kenya") %>%
  right_join(select(kenya_surv, -surv, -Lx)) %>%
  arrange(region, year, age) %>%
  mutate(fert = ifelse(is.na(fert), 0, fert))


# create constants
n_ages <- length(unique(kenya_fert$age))
n_periods <- length(unique(kenya_fert$year))
n_obs <- 1 # number of population observations (min. 1 for model to work)
# note: set t_obs = as.array(0) below if no actual observations

years <- unique(kenya_fert$year)
years <- c(years, max(years + 5))
ages <- c(unique(kenya_surv$age), Inf)

shona_mig <- reg %>%
  filter(basic_gender == 2) %>%
  filter(!is.na(year_migration_m) & !is.na(age_migration_m)) %>%
  mutate(
    mig_year = cut(year_migration_m, breaks = years, labels = head(years, -1), include.lowest = T, right = F),
    age = cut(basic_age, breaks = ages, labels = head(ages, -1), include.lowest = T, right = F)
  ) %>%
  count(mig_year, age) %>%
  mutate(mig_year = as.numeric(as.character(mig_year)), age = as.numeric(as.character(age))) %>% 
  complete(mig_year = head(years, -1), age = head(ages, -1), fill = list(n = 0))

pop_mig_obs <- matrix(
  data = shona_mig$n,
  nrow = n_periods,
  ncol = n_ages,
  byrow = TRUE
) %>% array(dim = c(1, dim = c(n_periods, n_ages)))

# shona_mig_mat <- matrix(shona_mig$n, nrow = n_periods, ncol = n_ages, byrow = TRUE)

shona_native <- reg %>%
  filter(basic_gender == 2) %>%
  filter(born_kenya_m == 1) %>%
  mutate(
    age = cut(basic_age, breaks = ages, labels = head(ages, -1), include.lowest = T, right = F)
  ) %>%
  count(age) %>%
  mutate(age = as.numeric(as.character(age))) %>%
  complete(age = head(ages, -1), fill = list(n = 0))

pop_native_obs <- matrix(
  data = shona_native$n,
  nrow = 1,
  ncol = n_ages
)

ZWE_props <- all_pop %>%
  filter(region == "Zimbabwe") %>%
  group_by(period) %>%
  mutate(pop_prop = pop/sum(pop)) %>%
  ungroup() %>%
  filter(period %in% head(years, -1), age %in% ages) %>%
  arrange(period, age)

mig_init_props <- matrix(
  ZWE_props$pop_prop,
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

# prepare population observation
shona_pop <- hhm %>%
  filter(gender == 2, !is.na(age)) %>%
  mutate(age = cut(as.numeric(age), breaks = c(seq(0, max(kenya_surv$age), by = 5), Inf), labels = seq(0, max(kenya_surv$age), by = 5), include_lowest = TRUE, right = FALSE)) %>%
  group_by(age) %>%
  summarise(n = n()) %>%
  # complete age groups (possibly missing from survey) 
  full_join(data.frame(age = as_factor(seq(0, max(kenya_surv$age), by = 5))), by = "age") %>%
  mutate(n = ifelse(is.na(n), 0, n))

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
  t_obs = as.array(15),
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

mod <- cmdstan_model(here("script/pop_reconstruct_disaggregated.stan"))


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

write_rds(to_save, here("output/intermediate/shona_disaggregated_model.rds"))

