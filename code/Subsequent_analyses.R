# Probit regression models predicting dropout from factor or mean scores
library(OpenMx)
library(dplyr)
library(here)

source(here("code/helper_functions.R"))

models <- c("AO-2S-PA", "MIM-2S-PA", "AO", "MIM", "Mean Score")
terms_2spa <- c("one_to_dropout_f", "latent_self_efficacy_to_dropout_f", "ELS_to_dropout_f")
terms_pa <- c("one_to_dropout_f", "self_efficacy_to_dropout_f", "ELS_to_dropout_f")

dir.create("rds", showWarnings = FALSE)
if (!file.exists("rds/score_df_continuous.rds")) {
  source("code/AO_MIM_Steps2_3.R")
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
df_ao <- cbind(df[, c("dropout_f", "ELS")], "self_efficacy" = df$self_efficacy_ao_c)
df_mim <- cbind(df[, c("dropout_f", "ELS")], "self_efficacy" = df$self_efficacy_mim_c)
df_mn <- cbind(df[, c("dropout_f", "ELS")], "self_efficacy" = df$self_efficacy_mn_c)

# 1 - NO RELIABILITY ADJUSTMENT (PA) ####

# using the centered FS variable from alignment optimization (AO) ####
mx_ao <- mxModel("mx_ao", 
  type = "RAM", 
  mxData(observed = df_ao, type = "raw"), 
  manifestVars = c("self_efficacy", "dropout_f", "ELS"), 
  # measurement error variance
  mxPath(from = "self_efficacy", arrows = 2, 
         free = FALSE, 
         labels = "data.self_efficacy_ev"), 
  # regression
  mxPath(from = "self_efficacy", to = "dropout_f", 
         values = 0, 
         labels = "self_efficacy_to_dropout_f"),  
  mxPath(from = "ELS", to = "dropout_f", 
         labels = "ELS_to_dropout_f"),
  # threshold for y
  mxThreshold("dropout_f", nThresh = 1),  
  # latent response scale
  mxPath(from = "dropout_f", arrows = 2, 
         free = FALSE, values = 1, 
         label ="dropout_f_with_dropout_f"), 
  # predictor variance - freely estimated
  mxPath(from = "self_efficacy", arrows = 2, 
         free = TRUE, values = 1, 
         labels = "self_efficacy_with_self_efficacy"),  
  mxPath(from = "ELS", arrows = 2, 
         free = TRUE, 
         labels = "ELS_with_ELS"),
  mxPath(from = "ELS", to = "self_efficacy", arrows = 2, 
         free = TRUE, values = 0,
         labels = "ELS_with_self_efficacy"), 
  mxPath(from = "one", to = "dropout_f", arrows = 1, 
         free = TRUE, values = 0, 
         label ="one_to_dropout_f"),
  mxPath(from = "one", to = "self_efficacy", 
         free = TRUE, values = 0,
         label = "one_to_self_efficacy"),
  mxPath(from = "one", to = "ELS", 
         label = "one_to_ELS"),
  # standardize the coefficient for self-efficacy
  mxAlgebra(self_efficacy_to_dropout_f * 
              sqrt(self_efficacy_with_self_efficacy) / 
              sqrt(dropout_f_with_dropout_f), name = "beta1"),
  mxCI("beta1"),
  mxFitFunctionML()
)
ao_fit <- mxTryHardOrdinal(mx_ao, intervals = TRUE)
ao_s <- summary(ao_fit)
ao_est_se_ci <- extract_est_se_ci(ao_fit, ao_s, terms_pa)

# using the centered FS variable from measurement invariance modeling (MIM) ####
mx_mim <- mxModel("mx_mim", 
                 type = "RAM", 
                 mxData(observed = df_mim, type = "raw"), 
                 manifestVars = c("self_efficacy", "dropout_f", "ELS"), 
                 # measurement error variance
                 mxPath(from = "self_efficacy", arrows = 2, 
                        free = FALSE, 
                        labels = "data.self_efficacy_ev"), 
                 # regression
                 mxPath(from = "self_efficacy", to = "dropout_f", 
                        values = 0, 
                        labels = "self_efficacy_to_dropout_f"),  
                 mxPath(from = "ELS", to = "dropout_f", 
                        labels = "ELS_to_dropout_f"),
                 # threshold for y
                 mxThreshold("dropout_f", nThresh = 1),  
                 # latent response scale
                 mxPath(from = "dropout_f", arrows = 2, 
                        free = FALSE, values = 1, 
                        label ="dropout_f_with_dropout_f"), 
                 # predictor variance - freely estimated
                 mxPath(from = "self_efficacy", arrows = 2, 
                        free = TRUE, values = 1, 
                        labels = "self_efficacy_with_self_efficacy"),  
                 mxPath(from = "ELS", arrows = 2, 
                        free = TRUE, 
                        labels = "ELS_with_ELS"),
                 mxPath(from = "ELS", to = "self_efficacy", arrows = 2, 
                        free = TRUE, values = 0,
                        labels = "ELS_with_self_efficacy"), 
                 mxPath(from = "one", to = "dropout_f", arrows = 1, 
                        free = TRUE, values = 0, 
                        label ="one_to_dropout_f"),
                 mxPath(from = "one", to = "self_efficacy", 
                        free = TRUE, values = 0,
                        label = "one_to_self_efficacy"),
                 mxPath(from = "one", to = "ELS", 
                        label = "one_to_ELS"),
                 # standardize the coefficient for self-efficacy
                 mxAlgebra(self_efficacy_to_dropout_f * 
                             sqrt(self_efficacy_with_self_efficacy) / 
                             sqrt(dropout_f_with_dropout_f), name = "beta1"),
                 mxCI("beta1"),
                 mxFitFunctionML()
)
mim_fit <- mxTryHardOrdinal(mx_mim, intervals = TRUE)
mim_s <- summary(mim_fit)
mim_est_se_ci <- extract_est_se_ci(mim_fit, mim_s, terms_pa)

# using the centered composite mean score variable (CMS) ####
mx_mn <- mxModel("mx_mn", 
                  type = "RAM", 
                  mxData(observed = df_mn, type = "raw"), 
                  manifestVars = c("self_efficacy", "dropout_f", "ELS"), 
                  # measurement error variance
                  mxPath(from = "self_efficacy", arrows = 2, 
                         free = FALSE, 
                         labels = "data.self_efficacy_ev"), 
                  # regression
                  mxPath(from = "self_efficacy", to = "dropout_f", 
                         values = 0, 
                         labels = "self_efficacy_to_dropout_f"),  
                  mxPath(from = "ELS", to = "dropout_f", 
                         labels = "ELS_to_dropout_f"),
                  # threshold for y
                  mxThreshold("dropout_f", nThresh = 1),  
                  # latent response scale
                  mxPath(from = "dropout_f", arrows = 2, 
                         free = FALSE, values = 1, 
                         label ="dropout_f_with_dropout_f"), 
                  # predictor variance - freely estimated
                  mxPath(from = "self_efficacy", arrows = 2, 
                         free = TRUE, values = 1, 
                         labels = "self_efficacy_with_self_efficacy"),  
                  mxPath(from = "ELS", arrows = 2, 
                         free = TRUE, 
                         labels = "ELS_with_ELS"),
                  mxPath(from = "ELS", to = "self_efficacy", arrows = 2, 
                         free = TRUE, values = 0,
                         labels = "ELS_with_self_efficacy"), 
                  mxPath(from = "one", to = "dropout_f", arrows = 1, 
                         free = TRUE, values = 0, 
                         label ="one_to_dropout_f"),
                  mxPath(from = "one", to = "self_efficacy", 
                         free = TRUE, values = 0,
                         label = "one_to_self_efficacy"),
                  mxPath(from = "one", to = "ELS", 
                         label = "one_to_ELS"),
                  # standardize the coefficient for self-efficacy
                  mxAlgebra(self_efficacy_to_dropout_f * 
                              sqrt(self_efficacy_with_self_efficacy) / 
                              sqrt(dropout_f_with_dropout_f), name = "beta1"),
                  mxCI("beta1"),
                  mxFitFunctionML()
)
mn_fit <- mxTryHardOrdinal(mx_mn, intervals = TRUE)
mn_s <- summary(mn_fit)
mn_est_se_ci <- extract_est_se_ci(mn_fit, mn_s, terms_pa)

# 2 - WITH RELIABILITY ADJUSTMENT (2S-PA) ####
# drop rows with missing reliability
df_a <- cbind(df[!is.na(df$approx_rel), ], 
              "self_efficacy" = df[!is.na(df$approx_rel),"self_efficacy_ao_c"],
              "self_efficacy_ev" = df[!is.na(df$approx_rel),"approx_ev"])
df_p <- cbind(df[!is.na(df$partial_rel), ], 
              "self_efficacy" = df[!is.na(df$partial_rel),"self_efficacy_mim_c"], 
              "self_efficacy_ev" = df[!is.na(df$partial_rel),"partial_ev"])


# using the centered FS variable from alignment optimization, with reliability correction (AO-2S-PA) ####
mx_tspa_ao <- mxModel("mx_tspa_ao", 
  type = "RAM", 
  mxData(observed = df_a, type = "raw"), 
  manifestVars = c("self_efficacy", "dropout_f", "ELS"), 
  latentVars = "latent_self_efficacy",
  # loading = reliability (as definition variable) = 1 for Bartlett FS
  mxPath(from = "latent_self_efficacy", to = "self_efficacy", 
         free = FALSE, values = 1, 
         labels = "latent_self_efficacy_to_self_efficacy"),
  # measurement error variance
  mxPath(from = "self_efficacy", arrows = 2, 
         free = FALSE, 
         labels = "data.self_efficacy_ev"), 
  # latent regression
  mxPath(from = "latent_self_efficacy", to = "dropout_f", 
         values = 0, 
         labels = "latent_self_efficacy_to_dropout_f"),  
  mxPath(from = "ELS", to = "dropout_f", 
         labels = "ELS_to_dropout_f"),
  # threshold for y
  mxThreshold("dropout_f", nThresh = 1),  
  # latent response scale
  mxPath(from = "dropout_f", arrows = 2, 
         free = FALSE, values = 1, 
         label ="dropout_f_with_dropout_f"), 
  # latent predictor variance - freely estimated
  mxPath(from = "latent_self_efficacy", arrows = 2, 
         free = TRUE, values = 1, 
         labels = "latent_self_efficacy_with_latent_self_efficacy"),  
  mxPath(from = "ELS", arrows = 2, 
         free = TRUE, 
         labels = "ELS_with_ELS"),
  mxPath(from = "ELS", to = "latent_self_efficacy", arrows = 2, 
         free = TRUE, values = 0,
         labels = "ELS_with_latent_self_efficacy"), 
  # latent intercept 
  mxPath(from = "one", to = c("latent_self_efficacy"), 
         free = FALSE, values = 0,
         label = "one_to_latent_self_efficacy"),  
  mxPath(from = "one", to = "dropout_f", arrows = 1, 
         free = TRUE, values = 0, 
         label ="one_to_dropout_f"),
  mxPath(from = "one", to = "self_efficacy", 
         free = TRUE, values = 0,
         label = "one_to_self_efficacy"),
  mxPath(from = "one", to = "ELS", 
         label = "one_to_ELS"),
  # standardize the coefficient for latent self-efficacy
  mxAlgebra(latent_self_efficacy_to_dropout_f * 
              sqrt(latent_self_efficacy_with_latent_self_efficacy) / 
              sqrt(dropout_f_with_dropout_f), name = "beta1"),
  mxCI("beta1"),
  mxFitFunctionML()
)

tspa_ao_fit <- mxTryHardOrdinal(mx_tspa_ao, intervals = TRUE)
tspa_ao_s <- summary(tspa_ao_fit)
tspa_ao_est_se_ci <- extract_est_se_ci(tspa_ao_fit, tspa_ao_s, terms_2spa)


# using the centered FS variable from measurement invariance modeling, with reliability correction (MIM-2S-PA) ####
mx_tspa_mim <- mxModel("mx_tspa_mim", 
  type = "RAM", 
  mxData(observed = df_p, type = "raw"), 
  manifestVars = c("self_efficacy", "dropout_f", "ELS"), 
  latentVars = "latent_self_efficacy",
  # loading = reliability (as definition variable) = 1 for Bartlett FS
  mxPath(from = "latent_self_efficacy", to = "self_efficacy", 
         free = FALSE, values = 1, 
         labels = "latent_self_efficacy_to_self_efficacy"),
  # measurement error variance
  mxPath(from = "self_efficacy", arrows = 2, 
         free = FALSE, 
         labels = "data.self_efficacy_ev"), 
  # latent regression
  mxPath(from = "latent_self_efficacy", to = "dropout_f", 
         values = 0, 
         labels = "latent_self_efficacy_to_dropout_f"),  
  mxPath(from = "ELS", to = "dropout_f", 
         labels = "ELS_to_dropout_f"),
  # threshold for y
  mxThreshold("dropout_f", nThresh = 1),  
  # latent response scale
  mxPath(from = "dropout_f", arrows = 2, 
         free = FALSE, values = 1, 
         label ="dropout_f_with_dropout_f"), 
  # latent predictor variance - freely estimated
  mxPath(from = "latent_self_efficacy", arrows = 2, 
         free = TRUE, values = 1, 
         labels = "latent_self_efficacy_with_latent_self_efficacy"),  
  mxPath(from = "ELS", arrows = 2, 
         free = TRUE, 
         labels = "ELS_with_ELS"),
  mxPath(from = "ELS", to = "latent_self_efficacy", arrows = 2, 
         free = TRUE, values = 0,
         labels = "ELS_with_latent_self_efficacy"), 
  # latent intercept 
  mxPath(from = "one", to = c("latent_self_efficacy"), 
         free = FALSE, values = 0,
         label = "one_to_latent_self_efficacy"),  
  mxPath(from = "one", to = "dropout_f", arrows = 1, 
         free = TRUE, values = 0, 
         label ="one_to_dropout_f"),
  mxPath(from = "one", to = "self_efficacy", 
         free = TRUE, values = 0,
         label = "one_to_self_efficacy"),
  mxPath(from = "one", to = "ELS", 
         label = "one_to_ELS"),
  # standardize the coefficient for latent self-efficacy
  mxAlgebra(latent_self_efficacy_to_dropout_f * 
              sqrt(latent_self_efficacy_with_latent_self_efficacy) / 
              sqrt(dropout_f_with_dropout_f), name = "beta1"),
  mxCI("beta1"),
  mxFitFunctionML()
)

tspa_mim_fit <- mxTryHardOrdinal(mx_tspa_mim, intervals = TRUE)
tspa_mim_s <- summary(tspa_ao_fit)
tspa_mim_est_se_ci <- extract_est_se_ci(tspa_mim_fit, tspa_mim_s, terms_2spa)

# Examine estimates ####
probit_table <- cbind(format_est_se_ci(ao_est_se_ci, 3), 
                      format_est_se_ci(mim_est_se_ci, 3),
                      format_est_se_ci(mn_est_se_ci, 3),
                      format_est_se_ci(tspa_ao_est_se_ci, 3), 
                      format_est_se_ci(tspa_mim_est_se_ci, 3))

ao_f_ind <-  c(ao_s$AIC, ao_s$BIC, ao_s$Minus2LogLikelihood) 
mim_f_ind <- c(mim_s$AIC, mim_s$BIC, mim_s$Minus2LogLikelihood) 
mn_f_ind <- c(mn_s$AIC, mn_s$BIC, mn_s$Minus2LogLikelihood) 
tspa_ao_f_ind <- c(tspa_ao_s$AIC, tspa_ao_s$BIC, tspa_ao_s$Minus2LogLikelihood) 
tspa_mim_f_ind <- c(tspa_mim_s$AIC, tspa_mim_s$BIC, tspa_mim_s$Minus2LogLikelihood) 

fit_ind <- rnd(cbind(ao_f_ind, mim_f_ind, mn_f_ind, tspa_ao_f_ind, tspa_mim_f_ind), 3)

probit_table <- rbind(probit_table,
                      c(ao_s$numObs, mim_s$numObs, mn_s$numObs, 
                        tspa_ao_s$numObs, tspa_mim_s$numObs),
                      fit_ind)
rownames(probit_table) <- c(rep(c("Est", "SE", "CI"), times = 3), "N", "AIC", "BIC", "-2LL")
colnames(probit_table) <- c(models[c(3:4)], "CMS",  models[c(1, 2)])
saveRDS(probit_table, 'rds/probit_reg_table_std.rds')






