# HERV_Bayesian_analyses
Scripts and data for Bayesian analyses in HERV paper ("Human endogenous retroviruses form a reservoir of T cell targets in hematological cancers with low mutational burden")

### Data
Original data are in excel files in the folder `data/raw`. These give the number of T-cell populations recognizing HERV or viral peptides for different patients, for individual HLA alleles present in those individuals.

### Scripts for analysis
The folder `scripts/R`contains one script for pre-processing the raw data (extracting relevant information, creating R-objects that can be used for the statistical analysis), and one script for running the [Stan](https://mc-stan.org) analyses via [RStan](https://mc-stan.org/users/interfaces/rstan.html). All Stan models that are called by this R-script are in the folder `scripts/Stan`. 

### Output
Processed data will be written to .rds files in the folder `data/processed`. The `stanfit` objects resulting from running the Stan models will be written to .rds files in the folder `results/Stan_fit`.
