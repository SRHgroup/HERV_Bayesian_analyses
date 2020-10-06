# HERV_Bayesian_analyses
Scripts and data for Bayesian analyses in HERV paper ("Human endogenous retroviruses form a reservoir of T cell targets in hematological cancers with low mutational burden")

### Data
Original data are in excel files in the folder `data/raw`. These give the number of T-cell populations recognizing HERV or viral peptides for different individuals (patients or controls), for the different HLA alleles present in those individuals.

### Scripts for analysis
The folder `scripts/R` contains three consecutively numbered scripts that will redo the analysis. They should be run in numerical order:
* **01_data_processing_github.R:** Script for pre-processing of raw data. Relevant information is extracted from the excel files, reformatted and saved as R-objects that can be used for the statistical analysis.
* **02_run_stan_models_github.R:** Script for running [Stan](https://mc-stan.org) analyses via [RStan](https://mc-stan.org/users/interfaces/rstan.html). All Stan scripts that are called by this R-script are in the folder `scripts/Stan`. 
* **03_plot_bayes_figures_github.R:** Script for plotting results.

### Output
Processed data will be written to .rds files in the folder `data/processed`. The `stanfit` objects resulting from running the Stan models will be written to .rds files in the folder `results/Stan_fit`. Figures created by the plotting script will be written to `results/figures`.
