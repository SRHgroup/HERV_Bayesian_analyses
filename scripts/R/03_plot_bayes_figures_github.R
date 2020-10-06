library(tidyverse)
library(tidybayes)

########################################################################################
# Function that creates list of ggplot annotate segments forming a significance bracket 
# for plotting. This list of annotate segments an be added to a ggplot object

bracket = function (x1, x2, y1, y2, lab, textpos = 0.5, endlen = 0.05, reverse = FALSE) {
    # textpos: fraction of diff from smallest x or y value
    # Vertical
    if (reverse) {
        endlen = - 1 * endlen
    }
    if (x1 == x2) {
        if (y2 < y1) {
            tmp = y1
            y1 = y2
            y2 = tmp
        }
        x1a = x1 - endlen; y1a = y1
        x2a = x2 - endlen; y2a = y2
        xtext = x1 - endlen; ytext = y1 + (y2 - y1) * textpos
        if (reverse) {
            vj = "center"; hj = "left"
        } else {
            vj = "center"; hj = "right"
        }
    # Horizontal
    } else if (y1 == y2) {
        if (x2 < x1) {
            tmp = x1
            x1 = x2
            x2 = tmp
        }
        x1a = x1; y1a = y1 - endlen
        x2a = x2; y2a = y2 - endlen
        xtext = x1 + (x2 - x1) * textpos; ytext = y1 + endlen
        if (reverse) {
            vj = "top"; hj = "center"
        } else {
            vj = "bottom"; hj = "center"
        }
    } else {
        print("Bracket line must be vertical or horizontal")
    }
    return(
        list(
            annotate("segment", x=x1, xend=x2, y=y1, yend=y2),
            annotate("segment", x=x1, xend=x1a, y=y1, yend=y1a),
            annotate("segment", x=x2, xend=x2a, y=y2, yend=y2a),
            annotate("text", x=xtext, y=ytext, label=lab, vjust=vj, hjust=hj)
        )
    )
}

########################################################################################

# Figure 3A
# Plot posterior for simple model of positive responders

fit = readRDS("../../results/Stan_fit/fit_07.rds")

df = fit %>% 
    gather_draws(theta[group]) %>% 
    mutate("Class" = case_when(
        group == 1 ~ "control",
        group == 2 ~ "patient",
        group == 3 ~ "patientAZA"
    ))

# Compute posterior probabilities that will be rendered on plot
m = fit %>%
    as.matrix() %>%
    as.data.frame()
ntot = nrow(m)

# Probability that proportion in patients > prop in controls:
n = m %>% filter(`theta[2]`>`theta[1]`) %>% nrow
pt_gt_con = sprintf("%.1f%%", 100 * n / ntot)

# Probability that proportion in AZA > prop in controls:
n = m %>% filter(`theta[3]`>`theta[1]`) %>% nrow
ptaza_gt_con = sprintf("%.1f%%", 100 * n / ntot)

maxval = df$.value %>%
    max() %>%
    `*`(1.01)
ggplot(df, aes(x=.value, y=Class, fill=Class)) + 
    stat_eye(.width=c(0.5, 0.9), show.legend=FALSE)+
    labs(y="") +
    theme_classic() +
    theme(panel.border = element_rect(fill=NA, size=1.5),
          axis.text = element_text(size=12)) + 
    scale_x_continuous(breaks=seq(0, 0.7, 0.1), limits = c(0, maxval)) +
    scale_y_discrete(limits=c("patientAZA", "patient", "control"),
                     labels=c("Post-AZA", "Pre-AZA", "Healthy donors")) + 
    scale_fill_manual(values=c("grey65", "#E3CCE5", "#D87BB8")) +
    labs(x="Proportion with immune response",
         title="Proportion of individuals with immune response to ERV") +
    bracket(0.67, 0.67, 2, 3, pt_gt_con, endlen=0.01) + 
    bracket(0.685, 0.685, 1, 3, ptaza_gt_con, endlen=0.01, textpos=0.25) 

ggsave("../../results/figures/3A_proportions_people_ERV_eye.pdf", width=17, height=12, units="cm")

########################################################################################

# Figure 3C
# Posterior for ERV proportion regression model 
# (proportion of peptides being positive, corrected for alleles, not normalised)

fit = readRDS("../../results/Stan_fit/fit_05.rds")

df = fit %>% 
    gather_draws(`^rr_.*`, regex=TRUE) %>% 
    mutate(logval = log(.value))

# Compute table of estimates and posterior probabilities for plot
posterior_table = df %>% 
    ungroup() %>%
    group_by(.variable) %>% 
    summarise(median = median(logval), 
              CI_lower=quantile(logval, 0.05), 
              CI_upper=quantile(logval, 0.95), 
              P_gt_0 = sum(logval>0) / n())

# Create strings giving values of posterior probabilities
# Note: piping to sprintf uses trick with braces to avoid injecting twice
# https://stackoverflow.com/questions/42385010/using-the-pipe-and-dot-notation
pt_gt_con = posterior_table %>%
    filter(.variable == "rr_pt_con") %>%
    pull(P_gt_0) %>%
    { sprintf(fmt="%.1f%%", 100 * .) }  

ptaza_gt_con = posterior_table %>%
    filter(.variable == "rr_aza_con") %>%
    pull(P_gt_0) %>%
    { sprintf(fmt="%.1f%%", 100 * .) }  

ptaza_gt_pt = posterior_table %>%
    filter(.variable == "rr_aza_pt") %>%
    pull(P_gt_0) %>%
    { sprintf(fmt="%.1f%%", 100 * .) }  

# Find posterior medians for plotting brackets
median_pt_gt_con = posterior_table %>%
    filter(.variable == "rr_pt_con") %>%
    pull(median)

median_ptaza_gt_con = posterior_table %>%
    filter(.variable == "rr_aza_con") %>%
    pull(median)

median_ptaza_gt_pt = posterior_table %>%
    filter(.variable == "rr_aza_pt") %>%
    pull(median)

ggplot(df, aes(x=logval, y=.variable, fill=.variable)) + 
    stat_eye(.width=c(0.5, 0.9), show.legend=FALSE) + 
    geom_vline(xintercept = 0, lty=2, lwd=0.3, alpha=0.6, color="blue") + 
    labs(y="") +
    theme_classic() +
    theme(panel.border = element_rect(fill=NA, size=1.5),
          axis.text = element_text(size=12)) + 
    scale_y_discrete(limits=c("rr_aza_pt", "rr_aza_con", "rr_pt_con"),
                     labels=c("Post-AZA vs Pre-AZA", "Post-AZA vs Healthy", "Pre-AZA vs Healthy")) + 
    scale_fill_manual(values=c("#D87BB8", "#D87BB8", "#D87BB8")) +
    labs(x="Log fold change",
         title="Log fold change in positive ERV peptides",
         subtitle="Corrected for HLA allele") + 
    bracket(0, median_pt_gt_con, 3.4, 3.4, pt_gt_con, endlen=0.03) + 
    bracket(0, median_ptaza_gt_con, 2.4, 2.4, ptaza_gt_con, endlen=0.03) + 
    bracket(0, median_ptaza_gt_pt, 1.51, 1.51, ptaza_gt_pt, endlen=0.03, textpos = 0.6)

ggsave("../../results/figures/3C_ERV_logfold_change_regression.pdf", width=17, height=12, units="cm")

########################################################################################

# Figure 3D
# Proportion of people responding to viral peptides
fit = readRDS("../../results/Stan_fit/fit_07_viral.rds")

df = fit %>% 
    gather_draws(theta[group]) %>% 
    mutate("Class" = case_when(
        group == 1 ~ "control",
        group == 2 ~ "patient",
        group == 3 ~ "patientAZA"
    ))

# Compute posterior probabilities that will be rendered on plot
m = fit %>%
    as.matrix() %>%
    as.data.frame()
ntot = nrow(m)

# Probability that proportion in patients < prop in controls:
n = m %>% filter(`theta[2]`<`theta[1]`) %>% nrow
pt_lt_con = sprintf("%.1f%%", 100 * n / ntot)

# Probability that proportion in AZA < prop in controls:
n = m %>% filter(`theta[3]`<`theta[1]`) %>% nrow
ptaza_lt_con = sprintf("%.1f%%", 100 * n / ntot)

maxval = df$.value %>%
    max() %>%
    `*`(1.01)

ggplot(df, aes(x=.value, y=Class, fill=Class)) + 
    stat_eye(.width=c(0.5, 0.9), show.legend=FALSE)+
    labs(y="") +
    scale_x_continuous(breaks=seq(0, 1, 0.1), limits = c(0,maxval)) +
    theme_classic() +
    theme(panel.border = element_rect(fill=NA, size=1.5),
          axis.text = element_text(size=12)) + 
    scale_y_discrete(limits=c("patientAZA", "patient", "control"),
                     labels=c("Post-AZA", "Pre-AZA", "Healthy donors")) + 
    scale_fill_manual(values=c("grey65", "lightskyblue2", "royalblue2")) +
    labs(x="Proportion with immune response",
         title="Proportion of individuals with immune response to viral peptides") +
    bracket(0.1, 0.1, 2, 3, pt_lt_con, endlen=0.01, reverse=TRUE) +
    bracket(0.085, 0.085, 1, 3, ptaza_lt_con, endlen=0.01, textpos=0.25, reverse=TRUE) 

ggsave("../../results/figures/3D_proportions_people_vir_eye.pdf", width=17, height=12, units="cm")

########################################################################################

# Figure 3F
# Log fold change of VIR peptides for 3 groups
fit = readRDS("../../results/Stan_fit/fit_08.rds")

df = fit %>% 
    gather_draws(`^rr_.*`, regex=TRUE) %>% 
    mutate(logval = log(.value))

# Compute table of estimates and posterior probabilities for plot
posterior_table = df %>% 
    ungroup() %>%
    group_by(.variable) %>% 
    summarise(median = median(logval), 
              CI_lower=quantile(logval, 0.05), 
              CI_upper=quantile(logval, 0.95), 
              P_gt_0 = sum(logval>0) / n())

# Note: piping to sprintf uses trick with braces to avoid injecting twice
# https://stackoverflow.com/questions/42385010/using-the-pipe-and-dot-notation
pt_lt_con = posterior_table %>%
    filter(.variable == "rr_pt_con") %>%
    pull(P_gt_0) %>%
    { sprintf(fmt="%.1f%%", 100 * (1 - .)) }  

ptaza_lt_con = posterior_table %>%
    filter(.variable == "rr_aza_con") %>%
    pull(P_gt_0) %>%
    { sprintf(fmt="%.1f%%", 100 * (1 - .)) }  

ptaza_gt_pt = posterior_table %>%
    filter(.variable == "rr_aza_pt") %>%
    pull(P_gt_0) %>%
    { sprintf(fmt="%.1f%%", 100 * .) }  

# Find posterior medians for plotting brackets
median_pt_lt_con = posterior_table %>%
    filter(.variable == "rr_pt_con") %>%
    pull(median)

median_ptaza_lt_con = posterior_table %>%
    filter(.variable == "rr_aza_con") %>%
    pull(median)

median_ptaza_gt_pt = posterior_table %>%
    filter(.variable == "rr_aza_pt") %>%
    pull(median)

ggplot(df, aes(x=logval, y=.variable, fill=.variable)) + 
    stat_eye(.width=c(0.5, 0.9), show.legend=FALSE) + 
    geom_vline(xintercept = 0, lty=2, lwd=0.3, alpha=0.6, color="blue") + 
    labs(y="") +
    theme_classic() +
    theme(panel.border = element_rect(fill=NA, size=1.5),
          axis.text = element_text(size=12)) + 
    scale_y_discrete(limits=c("rr_aza_pt", "rr_aza_con", "rr_pt_con"),
                     labels=c("Post-AZA vs Pre-AZA", "Post-AZA vs Healthy", "Pre-AZA vs Healthy")) + 
    scale_fill_manual(values=c("royalblue2", "royalblue2", "royalblue2")) +
    labs(x="Log fold change",
         title="Log fold change in positive viral peptides") +
    bracket(median_pt_lt_con, 0, 3.46, 3.46, pt_lt_con, endlen=0.03) +
    bracket(median_ptaza_lt_con, 0, 2.52, 2.52, ptaza_lt_con, endlen=0.03) +
    bracket(0, median_ptaza_gt_pt, 1.45, 1.45, ptaza_gt_pt, endlen=0.03)

ggsave("../../results/figures/3F_viral_logfold_change.pdf", width=17, height=12, units="cm")

########################################################################################

# Figure 3G
# Compute ratio between ERV proportions and viral proportions for three classes
# Done by random sampling from two separate posteriors, and then computing/plotting
# from random combination 
fitvir = readRDS("../../results/Stan_fit/fit_08.rds")
fiterv = readRDS("../../results/Stan_fit/fit_05.rds")

dfvir = fitvir %>% 
    spread_draws(`^theta.*`, regex=TRUE) %>%
    rename(p_vir_con = `theta[1]`, p_vir_pt = `theta[2]`, p_vir_aza = `theta[3]`)

dferv = fiterv %>% 
    spread_draws(`^p_.*`, regex=TRUE) %>%
    rename(p_erv_con = p_control, p_erv_pt = p_patient, p_erv_aza = p_aza)

dfcomb = bind_cols(dferv, dfvir) %>%
    mutate(
        logfc_pt_con_norm = log( (p_erv_pt / p_vir_pt) / (p_erv_con / p_vir_con)),
        logfc_aza_con_norm = log( (p_erv_aza / p_vir_aza) / (p_erv_con / p_vir_con)),
        logfc_aza_pt_norm = log( (p_erv_aza / p_vir_aza) / (p_erv_pt / p_vir_pt))
    )

df = dfcomb %>%
    pivot_longer(cols=starts_with("logfc_")) %>%
    dplyr::select(name, value)

# Compute table of estimates and posterior probabilities for plot
posterior_table = df %>% 
    ungroup() %>%
    group_by(name) %>% 
    summarise(median = median(value), 
              CI_lower=quantile(value, 0.05), 
              CI_upper=quantile(value, 0.95), 
              P_gt_0 = sum(value>0) / n())

# Note: piping to sprintf uses trick with braces to avoid injecting twice
# https://stackoverflow.com/questions/42385010/using-the-pipe-and-dot-notation
pt_gt_con = posterior_table %>%
    filter(name == "logfc_pt_con_norm") %>%
    pull(P_gt_0) %>%
    { sprintf(fmt="%.1f%%", 100 * .) }  

ptaza_gt_con = posterior_table %>%
    filter(name == "logfc_aza_con_norm") %>%
    pull(P_gt_0) %>%
    { sprintf(fmt="%.1f%%", 100 * .) }  

ptaza_lt_pt = posterior_table %>%
    filter(name == "logfc_aza_pt_norm") %>%
    pull(P_gt_0) %>%
    { sprintf(fmt="%.1f%%", 100 * (1 - .)) }  

# Find posterior medians for plotting brackets
median_pt_gt_con = posterior_table %>%
    filter(name == "logfc_pt_con_norm") %>%
    pull(median)

median_ptaza_gt_con = posterior_table %>%
    filter(name == "logfc_aza_con_norm") %>%
    pull(median)

median_ptaza_gt_pt = posterior_table %>%
    filter(name == "logfc_aza_pt_norm") %>%
    pull(median)

ggplot(df, aes(x=value, y=name, fill=name)) + 
    stat_eye(.width=c(0.5, 0.9), show.legend=FALSE, ) + 
    geom_vline(xintercept = 0, lty=2, lwd=0.3, alpha=0.6, color="blue") + 
    labs(y="") +
    theme_classic() +
    theme(panel.border = element_rect(fill=NA, size=1.5),
          axis.text = element_text(size=12)) + 
    scale_y_discrete(limits=c("logfc_aza_pt_norm", "logfc_aza_con_norm", "logfc_pt_con_norm"),
                     labels=c("Post-AZA vs Pre-AZA", "Post-AZA vs Healthy", "Pre-AZA vs Healthy")) + 
    scale_fill_manual(values=c("grey65", "grey65", "grey65")) +
    labs(x="Log fold change, normalised",
         title="Log fold change in positive ERV peptides",
         subtitle = "Corrected for HLA allele, normalised to viral response") +
    bracket(0, median_pt_gt_con, 3.43, 3.43, pt_gt_con, endlen=0.03) + 
    bracket(0, median_ptaza_gt_con, 2.43, 2.43, ptaza_gt_con, endlen=0.03) +
    bracket(median_ptaza_gt_pt, 0, 1.5, 1.5, ptaza_lt_pt, endlen=0.03)

ggsave("../../results/figures/3G_normalised_logfold_change.pdf", width=17, height=12, units="cm")

########################################################################################


