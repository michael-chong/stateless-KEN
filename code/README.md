# Code for population reconstruction

This folder contains the R code and Stan model files used to fit population reconstruction models and produce estimates. 

The `shona` and `pemba` folders contain scripts to fit models to the respective survey data. The Stan files specifying the models (which can be used for either dataset) are contained in the `models` subfolder.

The models included are variants of a Bayesian cohort component population projection approach. Below are brief descriptions of each model. For additional details, please consult the summary report "Exploration of methods to estimate stateless populations".

The first three models, listed below, test different ways of specifying the initial population. 

* `pop_reconstruction_uniform_init.stan`
* `pop_reconstruction_dirichlet_init.stan`
* `pop_reconstruction_lognormal_init.stan`

These models treat the survey observation as independent, age-specific counts. They assume a single initial population, with no further migration throughout the period of interest. A prior is placed on the initial population size, while the age proportions in the initial population are given either a uniform prior, or a user-specified prior which can be parameterized using either a Dirichlet distribution or relative proportions with log-Normal distributions. 

Next, `pop_reconstruction_lognormal_init_lognormal_obs.stan` allows for incorporation of the survey observation as proportions and total counts separately.

Finally, `pop_reconstruction_disaggregated.stan` models the population according to their reported period of migration, which in principle should allow for more accurate reconstruction of the past population by accounting for migration throughout the period of interest.