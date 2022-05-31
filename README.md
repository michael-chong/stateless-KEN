# Stateless populations estimation from Kenyan surveys

This repository contains reconstructions of stateless populations in Kenya based on survey data from the UNHCR from 2016 regarding the Kenyan Pemba population and from 2019 regarding the Kenyan Shona population.

# Contents

The repository is organized into 3 main folders: `code/`, `data/`, and `output/`.

The `data` folder contains some of the data files necessary to reproduce the simulations. The models require some the survey observations and historical demographic rates, which we take to be those from the 2019 World Population Prospects. Currently, the files containing survey observations are not included in this repository. 

The `code` folder contains the R code and Stan model files to fit the model and produce estimates of the historical populations. For illustrative and comparison purposes, different model variations included, which reflect different assumptions about observations and population processes. Some brief details describing the models are given in the file `code/README.md`. For additional context, please consult Section 2 of the report "Exploration of methods to estimate stateless populations".

The `output` folder contains select results and intermediate files. Typically, the R scripts in the R folder will write intermediate files to `output/intermediate/`, but these are currently not included in the repository because they contain copies of the survey data. Plots in the `output/pemba/` and `output/shona/` subfolders are generated from the notebooks `code/pemba_plots.Rmd` and `code/shona_plots.Rmd` respectively. 

# Contact

Michael Chong (myc.chong@mail.utoronto.ca)
