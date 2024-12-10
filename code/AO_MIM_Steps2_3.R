
# Determine an approximate invariance model via AO or a partial invariance model
# via MIM (Step 2) and compute factor scores and standard errors (Step 3)

library(lavaan)
library(sirt)
library(dplyr)
library(kableExtra)

# if data file does not already exist, download and prepare data
dir.create("rds", showWarnings = FALSE)
if (!file.exists("rds/dat.rds")) {
  source("code/download_data.R")
  source("code/prepare_data.R")
} else {
  dat <- readRDS("rds/dat.rds")
}

source("code/helper_functions.R")

### Step 2: Alignment Optimization (AO) ####

# See 'AO_MIM_Step1_determine_configural_model.R' for the steps leading up to 
# the final configural model below.
cfa_config <-  '
  group: ELS
  math =~ NA * i1 + l2_1 * i2 + l3 * i3 + l4_1 * i4 + l5_1 * i5
  i1 ~ nu1_1 * 1
  i2 ~ nu2_1 * 1
  i3 ~ nu3 * 1
  i4 ~ nu4_1 * 1
  i5 ~ nu5_1 * 1
  i1 ~~ theta1_1 * i1
  i2 ~~ theta2_1 * i2
  i3 ~~ theta3 * i3
  i4 ~~ theta4_1 * i4
  i5 ~~ theta5_1 * i5
  i1 ~~ i2
  i2 ~~ cov3 * i3
  i2 ~~ i4
  math ~~ 1 * math
  math ~ 0 * 1      
     
  group: HSLS
  math =~ NA * i1 + l2_2 * i2 + l4_2 * i4 + l5_2 * i5
  i1 ~ nu1_2 * 1
  i2 ~ nu2_2 * 1
  i4 ~ nu4_2 * 1
  i5 ~ nu5_2 * 1
  i1 ~~ theta1_2 * i1
  i2 ~~ theta2_2 * i2
  i4 ~~ theta4_2 * i4
  i5 ~~ theta5_2 * i5
  i1 ~~ i2
  i2 ~~ i4
  math ~~ 1 * math
  math ~ 0 * 1  
'
if (!file.exists("rds/fit_config.rds")) {
  fit_config  <- cfa(cfa_config, data = dat, group = "sample", 
                     estimator = "MLR", missing = "FIML", se = "robust.mlr")
  saveRDS(fit_config, "rds/fit_config.rds")
} else {
  fit_config <- readRDS("rds/fit_config.rds")
}
est_config <- lavInspect(fit_config, what = "est")
config_lambda <- merge(est_config$ELS$lambda, est_config$HSLS$lambda,
                       by = "row.names", all = TRUE)
config_lambda_mat <- t(config_lambda[, -1])
config_nu <- merge(est_config$ELS$nu, est_config$HSLS$nu,
                   by = "row.names", all = TRUE)
config_nu_mat <- t(config_nu[, -1])

# We next obtain a weight matrix of the sample sizes for each item by sample.
m_items <- paste0("i", 1:5)
wgt_mat <- as.matrix(rbind(
  colSums(!is.na(dat[dat$sample == "ELS", m_items])),
  colSums(!is.na(dat[dat$sample == "HSLS", m_items]))
))

# We perform alignment using the `invariance.alignment` function from the R
# package `sirt` with the configural loading and intercept estimates.
aligned_par <- invariance.alignment(lambda = config_lambda_mat,
                                    nu = config_nu_mat,
                                    fixed = TRUE, # fix SD of first group to one
                                    wgt = sqrt(wgt_mat))
# We modify the configural model specification string to freely estimate the 
# latent mean and variance in both groups, and substitute in the aligned 
# estimates from the configural model.
cfa_align <- '
  group: ELS
  math =~ l1_1 * i1 + l2_1 * i2 + l3 * i3 + l4_1 * i4 + l5_1 * i5
  i1 ~ nu1_1 * 1
  i2 ~ nu2_1 * 1
  i3 ~ nu3 * 1
  i4 ~ nu4_1 * 1
  i5 ~ nu5_1 * 1

  i1 ~~ theta1_1 * i1
  i2 ~~ theta2_1 * i2
  i3 ~~ theta3 * i3
  i4 ~~ theta4_1 * i4
  i5 ~~ theta5_1 * i5
  i1 ~~ i2
  i2 ~~ cov3 * i3
  i2 ~~ i4

  math ~~ NA * math
  math ~ NA * 1

  group: HSLS
  math =~ l1_2 * i1 + l2_2 * i2 + l4_2 * i4 + l5_2 * i5
  i1 ~ nu1_2 * 1
  i2 ~ nu2_2 * 1
  i4 ~ nu4_2 * 1
  i5 ~ nu5_2 * 1

  i1 ~~ theta1_2 * i1
  i2 ~~ theta2_2 * i2
  i4 ~~ theta4_2 * i4
  i5 ~~ theta5_2 * i5
  i1 ~~ i2
  i2 ~~ i4

  math ~~ NA * math
  math ~ NA * 1
'
cfa_align <- sub("l1_1", paste0(aligned_par$lambda.aligned[,1][1],
                                collapse = ","), cfa_align)
cfa_align <- sub("l1_2", paste0(aligned_par$lambda.aligned[,1][2],
                                collapse = ","), cfa_align)
cfa_align <- gsub("nu1_1", paste0(aligned_par$nu.aligned[,1][1],
                                  collapse = ","), cfa_align)
cfa_align <- gsub("nu1_2", paste0(aligned_par$nu.aligned[,1][2],
                                  collapse = ","), cfa_align)
fit_align  <- cfa(model = cfa_align,
                  data = dat,  estimator = "MLR", group = "sample",
                  missing = "FIML", se = "robust.mlr")
# obtain parameter estimates
est_align <- lavInspect(fit_align, what = "est")


### Step 2: Measurement Invariance Modeling (MIM) ####

# See 'MIM_Step2_detailed.R' for the MIM procedure resulting in the final 
# partial invariance model below
cfa_partial <-   '
  group: ELS
  math =~ NA * i1 + i2 + l3 * i3 + l4 * i4 + l5 * i5
  i1 ~ NA * 1
  i2 ~ NA * 1
  i3 ~ nu3 * 1
  i4 ~ nu4 * 1
  i5 ~ nu5 * 1

  i1 ~~ theta1_1 * i1
  i2 ~~ theta2_1 * i2
  i3 ~~ theta3 * i3
  i4 ~~ theta4_1 * i4
  i5 ~~ theta5_1 * i5
  i1 ~~ i2
  i2 ~~ cov3 * i3
  i2 ~~ i4

  math ~~ 1 * math
  math ~ 0 * 1

  group: HSLS
  math =~ NA * i1 + i2 + l4 * i4 + l5 * i5
  i1 ~ NA * 1
  i2 ~ NA * 1
  i4 ~ nu4 * 1
  i5 ~ nu5 * 1
  
  i1 ~~ theta1_2 * i1
  i2 ~~ theta2_2 * i2
  i4 ~~ theta4_2 * i4
  i5 ~~ theta5_2 * i5
  i1 ~~ i2
  i2 ~~ i4

  math ~~ NA * math
  math ~ NA * 1
'
fit_partial  <- cfa(cfa_partial, data = dat, group = "sample",
                    estimator = "MLR", missing = "FIML", se = "robust.mlr")
# extract cfa parameter estimates
est_partial <- lavInspect(fit_partial, what = "est")




### Step 3: Computation of factor scores (FS) and FS standard errors (SE) ####
fs_partial <- lavPredict(fit_partial, method = "bartlett", se = TRUE)
fs_align <- lavPredict(fit_align, method = "bartlett", se = TRUE)

# extract standard errors of the Bartlett scores
fs_partial_SE <- unlist(attributes(fs_partial)$se)
fs_align_SE <- unlist(attributes(fs_align)$se)

# store FS and FS SE in the score_df data frame
score_df <- as.data.frame(cbind(dat,
                                "partial_cont" = as.numeric(unlist(fs_partial)),
                                "partial_SE" = fs_partial_SE,
                                "approx_cont" = as.numeric(unlist(fs_align)),
                                "approx_SE" = fs_align_SE))

bartlett_rel <- function(psi, SE) { psi / (psi + SE^2) }

# compute reliability for each individual
rel_partial_ELS <- bartlett_rel(c(est_partial$ELS$psi),
                                score_df[score_df$sample == "ELS",
                                         "partial_SE"])
rel_partial_HSLS <- bartlett_rel(c(est_partial$HSLS$psi),
                                 score_df[score_df$sample == "HSLS",
                                          "partial_SE"])
score_df$partial_rel <- c(rel_partial_ELS, rel_partial_HSLS)

rel_approx_ELS <- bartlett_rel(c(est_align$ELS$psi),
                               score_df[score_df$sample == "ELS",
                                        "approx_SE"])
rel_approx_HSLS <- bartlett_rel(c(est_align$HSLS$psi),
                                score_df[score_df$sample == "HSLS",
                                         "approx_SE"])
score_df$approx_rel <- c(rel_approx_ELS, rel_approx_HSLS)

# store error variances
score_df <- as.data.frame(
  cbind(score_df, 
        "partial_ev" = score_df$partial_SE^2, 
        "approx_ev" = score_df$approx_SE^2))

score_df[, 3:11] <- apply(score_df[, 3:11], FUN = as.numeric, MARGIN = 2)

# create standardized versions of the score variables
score_df$self_efficacy_ao_s <- standardize_scores(score_df$approx_cont)
score_df$self_efficacy_mim_s <- standardize_scores(score_df$partial_cont)
score_df$self_efficacy_mn_s <- standardize_scores(score_df$mean_score)

saveRDS(score_df, "rds/score_df_continuous.rds")
saveRDS(est_partial, "rds/est_partial.rds")
saveRDS(est_align, "rds/est_align.rds")
saveRDS(fit_partial, "rds/fit_partial.rds")
saveRDS(fit_align, "rds/fit_align.rds")

# dMACS effect sizes
pinsearch::pin_effsize(fit_partial)
pinsearch::pin_effsize(fit_align)

# correlations
vars <- c("mean_score", "partial_cont", "approx_cont")
els_cor <- cor(score_df[score_df$sample == "ELS", vars], use = "complete.obs")
hsls_cor <- cor(score_df[score_df$sample == "HSLS", vars], use = "complete.obs")
all_cor <- cor(score_df[, vars], use = "complete.obs")
t_cor <- rbind(els_cor[upper.tri(els_cor)], 
               hsls_cor[upper.tri(els_cor)], 
               all_cor[upper.tri(els_cor)])
colnames(t_cor) <-c("partial-mean", "approx-mean", "approx-partial")
rownames(t_cor) <- c("ELS", "HSLS", "overall")
t_cor


## Build tables
build_mim_table(dig = 2)
build_ao_mim_parameter_table(dig = 2)
build_ao_mim_latent_table(dig = 2)
build_mnsd_table(dig = 2)

