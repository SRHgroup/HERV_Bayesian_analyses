# Script for running Stan models. Run after 01_data_processing_github.R
library(tidyverse)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

##############################################################

# Run stan model 05 (effect of cancer and aza accounting for alleles, sum to zero constraints)
df = readRDS("../../data/processed/all_erv_tidy.rds")

dat = list(
    N = nrow(df),
    N_group = df %>% 
        pull(group) %>% 
        n_distinct(),
    N_hla = df %>% 
        pull(allele_num) %>% 
        n_distinct(),
    group = df %>% 
        pull(group),
    hla = df %>% 
        pull(allele_num),
    n = df %>% 
        pull(n_tested),
    y = df %>% 
        pull(n_positive)
)

sm = stan_model("../Stan/model_05.stan")
fit = sampling(sm, data=dat, chains=3, iter=10000)
saveRDS(fit, "../../results/Stan_fit/fit_05.rds")

# Example of post-processing to compute probabilities of interest:
m = fit %>%
    as.matrix() %>%
    as.data.frame()

ntot = nrow(m)

# Probability that proportion in patients before AZA > prop in controls:
n = m %>% filter(p_patient>p_control) %>% nrow
n / ntot

# Probability that proportion in patients after AZA > prop in controls:
n = m %>% filter(p_aza>p_control) %>% nrow
n / ntot

# Probability that proportion after AZA > prop before AZA:
n = m %>% filter(p_aza>p_patient) %>% nrow
n / ntot

#########################################################################

# Run stan model for comparing grand total (not stratified on alleles)
# Estimating proportion of PEOPLE that respond to erv in each group

df = readRDS("../../data/processed/overall.rds")

dat = list(
    N = nrow(df),
    group = df %>% 
        pull(group),
    n = df %>% 
        pull(n_tested),
    y = df %>% 
        pull(n_positive)
)

sm = stan_model("../Stan/model_07.stan")
fit = sampling(sm, data=dat, chains=3, iter=10000)
saveRDS(fit, "../../results/Stan_fit/fit_07.rds")

#########################################################################
# Estimating proportion of people that respond to VIRAL peptides in each group

df = readRDS("../../data/processed/all_vir_tidy.rds") 

df2 = df %>% 
    group_by(group_num) %>% 
    summarise(n_tested = n(), n_pos = sum(n_positive>0))

dat = list(
    N = nrow(df2),
    group = df2 %>% 
        pull(group_num),
    n = df2 %>% 
        pull(n_tested),
    y = df2 %>% 
        pull(n_pos)
)

sm = stan_model("../Stan/model_07.stan")
fit = sampling(sm, data=dat, chains=3, iter=10000)
saveRDS(fit, "../../results/Stan_fit/fit_07_viral.rds")

#########################################################################
# Estimate proportion of viral peptides recognised for each grouping

df = readRDS("../../data/processed/all_vir_tidy.rds") 

dat = list(
    N = nrow(df),
    N_group = df %>%
        pull(group_num) %>%
        n_distinct(),
    group = df %>% 
        pull(group_num),
    n = df %>% 
        pull(n_tested),
    y = df %>% 
        pull(n_positive)
)

sm = stan_model("../Stan/model_08.stan")
fit = sampling(sm, data=dat, chains=3, iter=10000)
saveRDS(fit, "../../results/Stan_fit/fit_08.rds")

#########################################################################
# Logistic regression: HERV post aza standardized counts as predictors for clinical outcome
dat = readRDS("../../data/processed/dat_herv_post.rds")
sm2 = stan_model("../Stan/logreg_02.stan")
fit = sampling(sm2, data=dat, chains=3, iter=10000, control = list(adapt_delta = 0.95))
saveRDS(fit, "../../results/Stan_fit/fit_herv_post_logistic.rds")

#########################################################################
# Logistic regression: connection between VIR response and clinical outcome
dat = readRDS("../../data/processed/dat_vir_post.rds")
sm2 = stan_model("../Stan/logreg_02.stan")
fit = sampling(sm2, data=dat, chains=3, iter=10000, control = list(adapt_delta = 0.99))
saveRDS(fit, "../../results/Stan_fit/fit_vir_post_logistic.rds")

#########################################################################
# Logistic regression to predict clinical outcome
# Predictors are HERV, VIR, and interaction
# Note: binary predictors (not counts): presence or absence of response in patient
dat = readRDS("../../data/processed/dat_herv_vir_interaction.rds")
sm5 = stan_model("../Stan/logreg_05.stan")
fit = sampling(sm5, data=dat, chains=3, iter=10000, control = list(adapt_delta = 0.9))
saveRDS(fit, "../../results/Stan_fit/fit_interact_logistic.rds")




