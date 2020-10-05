# Script for extracting, reformatting, and saving relevant information from original excel files
# Resulting data is meant to be used in Stan scripts
# Somewhat messy, and could surely be simplified. 

library(tidyverse)
library(readxl)
library(tidybayes)

##########################################################################################
# Functions for processing data 
##########################################################################################

# Reformat raw dataframe read from excel file:
#    Add id column (person id)
#    Add class column (patient or control)
#    Convert to long format
#    Clean allele names (remove "HLA-") from new allele column (pivoted)
reformat_df = function (df, id_range, clss, valto) {
    df2 = df %>%
        mutate(id = id_range[1]:id_range[2], class=clss) %>%
        mutate(class_num = 
                   case_when(
                       class == "patient" ~ 1,
                       class == "control" ~ 0
                   )
        ) %>%
        pivot_longer(cols=starts_with("HLA"), names_to="allele", values_to=valto) %>%
        mutate(allele = str_replace(allele, "HLA-", "")) %>%
        mutate(allele_num = 
                   case_when(
                       allele == "A0101" ~ 1,
                       allele == "A0201" ~ 2,
                       allele == "B0702" ~ 3,
                       allele == "B0801" ~ 4
                   )
        )
    return(df2)
}

##########################################################################################
# Join reformatted data frames containing info on n_tested and n_positive, clean up
# Infer meaning of missing values (if testing done, then missing=0, else missing=NA)
join_df = function(df_ntest, df_npos) {
    df_res = left_join(df_ntest, df_npos,  
                       by=c("id", "allele", "allele_num", "class", "class_num")) %>%
        dplyr::select(id, class, class_num, allele, allele_num, n_tested, n_positive) %>%
        filter(!is.na(n_tested)) %>%
        mutate(n_positive = ifelse(is.na(n_positive), 0, n_positive))
}

##########################################################################################
# Read data from excel file, construct matrix with counts of recognised peptides for individual alleles
# Note: option standardize subtracts mean and divides by std.
#       This is used on predictors for logistic regression analysis
pep2mat = function(filename, range, standardize=FALSE) {
    mat = read_excel(filename, range=range, na="NA") %>%
        replace_na(list(A1=0,A2=0,B7=0,B8=0)) %>%
        as.matrix()
    if(standardize) {
        mat = mat %>% scale()
    }
    return(mat)
}


##########################################################################################
# Read data from excel file, construct binary vector indicating presence/absence of response
pep2vecbin = function(filename, range, standardize=FALSE) {
    vecbin = read_excel(filename, range=range, na="NA") %>%
        replace_na(list(A1=0,A2=0,B7=0,B8=0)) %>%
        mutate(pepsum = A1 + A2 + B7 + B8) %>%
        mutate(response = if_else(pepsum> 0, 1, 0)) %>%
        pull(response)
    if(standardize) {
        vecbin = vecbin %>% scale()
    }
    return(vecbin)
}


##########################################################################################
# Data processing 
##########################################################################################

# Create file "overall.rds" for models looking at grand total in groups
overall = read_excel("../../data/raw/HERV T cell data_04062020_v2.xlsx", sheet=2, range="A69:C72")
overall = overall %>%
    mutate(n_tested = Postive + Negative) %>%
    rename(n_positive = Postive) %>%
    mutate(group = case_when(
        Cohort == "Control" ~ 1,
        Cohort == "Before treatment" ~ 2,
        Cohort == "After treatment" ~ 3
    ))
saveRDS(overall, "../../data/processed/overall.rds")

###############################################################################################
# Process data for HERV. Create file "all_erv_tidy.rds"
# Join information about number tested and number positive
# Use to infer meaning of missing values (if testing done, then missing=0, else missing=NA)

pt_ntested = read_excel("../../data/raw/HERV T cell data_04062020_v2.xlsx", sheet=2, range="B2:F36")
pt_erv_before = read_excel("../../data/raw/HERV T cell data_04062020_v2.xlsx", sheet=2, range="G2:K36")
pt_erv_after = read_excel("../../data/raw/HERV T cell data_04062020_v2.xlsx", sheet=2, range="L2:P36")
con_ntested = read_excel("../../data/raw/HERV T cell data_04062020_v2.xlsx", sheet=2, range="B39:F66")
con_erv = read_excel("../../data/raw/HERV T cell data_04062020_v2.xlsx", sheet=2, range="G39:K66")

df1 = reformat_df(pt_ntested, c(1,34), "patient", "n_tested")
df2 = reformat_df(pt_erv_before, c(1,34), "patient", "n_positive")
df3 = reformat_df(pt_erv_after, c(1,34), "patient", "n_positive")
df4 = reformat_df(con_ntested, c(35,61), "control", "n_tested")
df5 = reformat_df(con_erv, c(35,61), "control", "n_positive")

pt_erv_before_tidy = join_df(df1, df2) %>% 
    filter(id != 2) %>%                                 # Patient 2: testing was not done
    mutate(aza=0)
pt_erv_after_tidy = join_df(df1, df3) %>%
    mutate(aza=1)
con_erv_tidy = join_df(df4, df5) %>%
    mutate(aza=0)

# Combine all data frames into one, save to disk
all_erv_tidy = bind_rows(pt_erv_before_tidy, pt_erv_after_tidy, con_erv_tidy) %>%
    dplyr::select(id, class, class_num, allele, allele_num, aza, n_tested, n_positive) %>%
    mutate(group = case_when(
        class == "control" ~ 1,
        class == "patient" & aza == 0 ~ 2,
        class == "patient" & aza == 1 ~ 3
        )
    )
saveRDS(all_erv_tidy, "../../data/processed/all_erv_tidy.rds")

###############################################################################################
# Process data for VIR (non-HERV, viral peptides) Create file "all_erv_tidy.rds"

vir_con_tot = read_excel("../../data/raw/HERV T cell data_04062020_v2.xlsx", sheet=2, 
                         range="L40:L66", col_names = "n_tested", na = "NA") %>% pull(n_tested)
vir_con_pos = read_excel("../../data/raw/HERV T cell data_04062020_v2.xlsx", sheet=2, 
                         range="M40:M66", col_names = "n_positive", na = "NA") %>% pull(n_positive)
vir_pt_tot = read_excel("../../data/raw/HERV T cell data_04062020_v2.xlsx", sheet=2, 
                        range="Q3:Q36", col_names = "n_tested", na = "NA") %>% pull(n_tested)
vir_pt_pos = read_excel("../../data/raw/HERV T cell data_04062020_v2.xlsx", sheet=2, 
                        range="R3:R36", col_names = "n_positive", na = "NA") %>% pull(n_positive)
vir_aza_pos = read_excel("../../data/raw/HERV T cell data_04062020_v2.xlsx", sheet=2, 
                         range="S3:S36", col_names = "n_positive", na = "NA") %>% pull(n_positive)

ptid = read_excel("../../data/raw/HERV T cell data_04062020_v2.xlsx", sheet=2, 
                  range="A3:A36", col_names = "id", na = "NA") %>% pull(id)
conid = read_excel("../../data/raw/HERV T cell data_04062020_v2.xlsx", sheet=2, 
                   range="A40:A66", col_names = "id", na = "NA") %>% pull(id)

vir_con = tibble(
    group = "Control",
    group_num = 1,
    id = conid,
    n_tested = vir_con_tot,
    n_positive = vir_con_pos
)

vir_pt = tibble(
    group = "Patient",
    group_num = 2,
    id = ptid,
    n_tested = vir_pt_tot,
    n_positive = vir_pt_pos
)

vir_aza = tibble(
    group = "Patient+AZA",
    group_num = 3,
    id = ptid,
    n_tested = vir_pt_tot,
    n_positive = vir_aza_pos
)

all_vir_tidy = bind_rows(vir_con, vir_pt, vir_aza) %>%
    drop_na()
saveRDS(all_vir_tidy, "../../data/processed/all_vir_tidy.rds")

###############################################################################################
# Prepare data for Stan analysis using logistic regression (responder ~ hla + pepcounts)
# Specifically Stan model will need the following:
# y: vector of outcome (responder / non-responder, coded as 1/0)
# hla: design matrix indicating which hla alleles individuals have
#      each row is an individual, the four columns contain indicators for hla being present
# herv_pre: design matrix for HERV response before AZA (number of peptides recognised)
#           0 can mean 0 or it can mean allele not present (taken care of by hla matrix)
# herv_post: same, after AZA
# vir_pre: as above but for viral antigens
# vir_pos: see above

y = read_excel("../../data/raw/Figure 4D-F_25082020.xlsx", range="B3:B36", col_names=c("responder")) %>%
    mutate(responder = if_else(responder == "Responder", 1, 0)) %>%
    pull(responder)

hla = read_excel("../../data/raw/Figure 4D-F_25082020.xlsx", range="C2:F36", col_names=TRUE) %>%
    mutate(across(everything(), ~if_else(. > 0, 1, 0, missing=0))) %>%
    as.matrix()

# connection between HERV response and clinical outcome
herv_post = pep2mat("../../data/raw/Figure 4D-F_25082020.xlsx", range="O2:R36", standardize=TRUE)
dat = list(N=length(y), y=y, hla=hla, pepcount=herv_post)
saveRDS(dat, "../../data/processed/dat_herv_post.rds")

# connection between VIR response and clinical outcome:
vir_post = pep2mat("../../data/raw/Figure 4D-F_25082020.xlsx", range="Y2:AB36", standardize=TRUE)
dat = list(N=length(y), y=y, hla=hla, pepcount=vir_post)
saveRDS(dat, "../../data/processed/dat_vir_post.rds")

# Interaction model: herv and vir are one column each (binary indicator of response or not)
# interaction is third predictor (herv * vir => 1 when both vir and herv response)
herv_post_bin = pep2vecbin("../../data/raw/Figure 4D-F_25082020.xlsx", range="O2:R36", standardize=FALSE)
vir_post_bin = pep2vecbin("../../data/raw/Figure 4D-F_25082020.xlsx", range="Y2:AB36", standardize=FALSE)
dat = list(N=length(y), y=y, hla=hla, herv=herv_post_bin, vir=vir_post_bin)
saveRDS(dat, "../../data/processed/dat_herv_vir_interaction.rds")



