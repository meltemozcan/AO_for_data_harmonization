# Probit regression models predicting dropout from standardized factor or mean scores
library(umx)
library(dplyr)
library(here)

source(here("code/helper_functions.R"))

models <- c("AO-2S-PA", "MIM-2S-PA", "AO", "MIM", "Mean Score")

dir.create("rds", showWarnings = FALSE)
if (!file.exists("rds/score_df_continuous.rds")) {
  source(here("code/AO_MIM_Steps1through3.R"))
  df <- readRDS(here("rds/score_df_continuous.rds"))
} else {
  df <- readRDS(here("rds/score_df_continuous.rds"))
}

# create a factor version of the dropout variable
df$dropout_f <- mxFactor(df$dropout, levels = c(0, 1))
# create a new variable to indicate sample membership
df$ELS <- 0
df$ELS[df$sample=="ELS"] <- 1 # ELS: 1, HSLS: 0

# prepare dfs with the standardized self_efficacy variables 
df_ao <- cbind(df[, c("dropout_f", "ELS")], "self_efficacy" = df$self_efficacy_ao_s)
df_mim <- cbind(df[, c("dropout_f", "ELS")], "self_efficacy" = df$self_efficacy_mim_s)
df_mn <- cbind(df[, c("dropout_f", "ELS")], "self_efficacy" = df$self_efficacy_mn_s)

# 1 - NO RELIABILITY ADJUSTMENT (PA) ####
umx_ao <- umxRAM("
    dropout_f ~ self_efficacy + ELS
    dropout_f ~ 0 * 1
    dropout_f ~~ 1 * dropout_f", 
    name = "probit regression (AO)", data = df_ao, tryHard = "yes")
umx_mim <- umxRAM(
   "dropout_f ~ self_efficacy + ELS
    dropout_f ~ 0 * 1
    dropout_f ~~ 1 * dropout_f", 
    name = "probit regression (MIM)", data = df_mim, tryHard = "yes")
umx_mn <- umxRAM(
   "dropout_f ~ self_efficacy + ELS
    dropout_f ~ 0 * 1
    dropout_f ~~ 1 * dropout_f", 
    name = "probit regression (Mean Scores)", data = df_mn, tryHard = "yes")

# equivalent model specification without using a model string:
umx_ao2 <- umxRAM("probit regression (AO)",
   data = df_ao, tryHard = "yes",
   umxPath(from = c("self_efficacy", "ELS"), to = "dropout_f"), # regression paths
   umxPath(from = "self_efficacy", to = "ELS", arrows = 2), # covariances
   umxPath(var = c("self_efficacy", "ELS")), # residual variances for predictors
   umxPath(means = c("self_efficacy", "ELS")), # means of predictors
   umxPath(v1m0 = "dropout_f") # variance of 1 and mean of zero
)
# Table: Parameter loadings for model 'probit regression'
#       name                             | Estimate|SE    |type          |
#   :--|:--------------------------------|--------:|:-----|:-------------|
#   5  |ELS_with_self_efficacy           |   -0.126|0.003 |Manifest Cov  |
#   1  |ELS_to_dropout_f                 |   -0.158|0.02  |Manifest path |
#   2  |self_efficacy_to_dropout_f       |   -0.163|0.011 |Manifest path |
#   7  |one_to_ELS                       |    0.408|0.002 |Mean          |
#   8  |one_to_self_efficacy             |   -0.017|0.006 |Mean          |
#   3  |dropout_f_with_dropout_f         |    1.000|0     |Residual      |
#   4  |ELS_with_ELS                     |    0.242|0.002 |Residual      |
#   6  |self_efficacy_with_self_efficacy |    1.002|0.008 |Residual      |
#   
#   Model Fit: χ²(1) = NA, p = NA; CFI = NA; TLI = NA; RMSEA = NA
#   Algebra'threshMat' = 1.364CI95[1.341, 1.387]. p-value < 0.001


# 2 - WITH RELIABILITY ADJUSTMENT (2S-PA) ####
# drop rows with missing reliability
df_a <- cbind(df[!is.na(df$approx_rel), ], 
              "self_efficacy" = df[!is.na(df$approx_rel),"self_efficacy_ao_s"],
              "self_efficacy_ev" = df[!is.na(df$approx_rel),"approx_ev"])
df_p <- cbind(df[!is.na(df$partial_rel), ], 
              "self_efficacy" = df[!is.na(df$partial_rel),"self_efficacy_mim_s"], 
              "self_efficacy_ev" = df[!is.na(df$partial_rel),"partial_ev"])

umx_tspa_ao <- umxRAM("2spa (AO)",
  umxPath(c("math", "ELS"), to = "dropout_f"), # main effects
  # loading = reliability (as definition variable) = 1 for Bartlett FS
  umxPath("math", to = "self_efficacy", values = 1, free = FALSE),
  # error variance (as definition variable), fs_SE^2
  umxPath(var = "self_efficacy", labels = "data.self_efficacy_ev", free = FALSE),
  umxPath(unique.pairs = c("math", "ELS")), # covariances of predictors
  umxPath(means = c("math", "ELS")),  # means of predictors
  umxPath(v1m0 = "dropout_f"), # for model identification, dropout_f ~ N(0, 1)
  data = df_a, tryHard = "yes"
)
umx_tspa_mim <- umxRAM("2spa (MIM)",
  umxPath(c("math", "ELS"), to = "dropout_f"), # main effects
  # loading = reliability (as definition variable) = 1 for Bartlett FS
  umxPath("math", to = "self_efficacy", values = 1, free = FALSE),
  # error variance (as definition variable), fs_SE^2
  umxPath(var = "self_efficacy", labels = "data.self_efficacy_ev", free = FALSE),
  umxPath(unique.pairs = c("math", "ELS")), # covariances of predictors
  umxPath(means = c("math", "ELS")), # means of predictors
  umxPath(v1m0 = "dropout_f"), # for model identification, dropout_f ~ N(0, 1)
  data = df_p, tryHard = "yes"
)



# Examine estimates ####

terms_2spa <- c("dropout_f_dev1", "math_to_dropout_f", "ELS_to_dropout_f")
terms_pa <- c("dropout_f_dev1", "self_efficacy_to_dropout_f", "ELS_to_dropout_f")

umx_est_ao <- get_path_ests_umx(terms_pa, umx_ao)
#                 Est     SE      LB      UB
# (intercept)   -1.3638 0.0118 -1.3871 -1.3407
# self-efficacy -0.1631 0.0114 -0.1855 -0.1408
# ELS           -0.1580 0.0199 -0.1971 -0.1189
umx_est_ao_tspa <- get_path_ests_umx(terms_2spa, umx_tspa_ao)
#                   Est     SE      LB      UB
# (intercept)   -1.4043 0.0134 -1.4307 -1.3781
# self-efficacy -0.1764 0.0123 -0.2006 -0.1525
# ELS           -0.1999 0.0234 -0.2460 -0.1541
umx_est_mim <- get_path_ests_umx(terms_pa, umx_mim)
#                   Est     SE      LB      UB
# (intercept)   -1.3636 0.0118 -1.3869 -1.3405
# self-efficacy -0.1617 0.0114 -0.1840 -0.1393
# ELS           -0.1578 0.0199 -0.1968 -0.1188
umx_est_mim_tspa <- get_path_ests_umx(terms_2spa, umx_tspa_mim)
#                   Est     SE      LB      UB
# (intercept)   -1.4040 0.0134 -1.4304 -1.3779
# self-efficacy -0.1745 0.0122 -0.1986 -0.1506
# ELS           -0.1997 0.0234 -0.2457 -0.1538
umx_est_mn <- get_path_ests_umx(terms_pa, umx_mn)
#                   Est     SE      LB      UB
# (intercept)   -1.3618 0.0118 -1.3851 -1.3386
# self-efficacy -0.1611 0.0114 -0.1837 -0.1386
# ELS           -0.1619 0.0200 -0.2013 -0.1226

probit_table <- cbind(f_ests_fit(umx_tspa_ao, terms_2spa, 3), 
                      f_ests_fit(umx_tspa_mim, terms_2spa, 3), 
                      f_ests_fit(umx_ao, terms_pa, 3),
                      f_ests_fit(umx_mim, terms_pa, 3), 
                      f_ests_fit(umx_mn, terms_pa, 3))
rownames(probit_table) <- c(rep(c("Est", "SE", "CI"), times = 3), "N", "AIC", "BIC", "-2LL")
colnames(probit_table) <- models
saveRDS(probit_table, 'rds/probit_reg_table_std.rds')

# Path diagrams ####
# plot(umx_mim)
# plot(umx_tspa_mim)
# plot(umx_ao)
# plot(umx_tspa_ao)
