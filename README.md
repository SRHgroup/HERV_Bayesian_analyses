# HERV_Bayesian_analyses
Scripts and data for Bayesian analyses in HERV paper:

[Human endogenous retroviruses form a reservoir of T cell targets in hematological cancers with low mutational burden](https://www.nature.com/articles/s41467-020-19464-8). Nature Communications, 2020, 11(1), 5660. https://doi.org/10.1038/s41467-020-19464-8

Saini, S. K., Ørskov, A. D., Bjerregaard, A.-M., Unnikrishnan, A., Holmberg-Thydén, S., Borch, A., Jensen, K. V., Anande, G., Bentzen, A. K., Marquard, A. M., Tamhane, T., Treppendahl, M. B., Gang, A. O., Dufva, I. H., Szallasi, Z., Ternette, N., Pedersen, A. G., Eklund, A. C., Pimanda, J., … Hadrup, S. R.


### Data
Original data are in excel files in the folder `data/raw`. These give the number of T-cell populations recognizing HERV or viral peptides for different individuals (patients or controls), for the different HLA alleles present in those individuals.

### Scripts for analysis
The folder `scripts/R` contains three consecutively numbered scripts that will redo the analysis. They should be run in numerical order:
* **01_data_processing_github.R:** Script for pre-processing of raw data. Relevant information is extracted from the excel files, reformatted and saved as R-objects that can be used for the statistical analysis.
* **02_run_stan_models_github.R:** Script for running [Stan](https://mc-stan.org) analyses via [RStan](https://mc-stan.org/users/interfaces/rstan.html). All Stan scripts that are called by this R-script are in the folder `scripts/Stan`. 
* **03_plot_bayes_figures_github.R:** Script for plotting results.

### Output
Processed data will be written to .rds files in the folder `data/processed`. The `stanfit` objects resulting from running the Stan models will be written to .rds files in the folder `results/Stan_fit`. Figures created by the plotting script will be written to `results/figures`.
