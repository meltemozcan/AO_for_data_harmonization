# Functions:

  # mnsd()
  # rnd()
  # f_ci()
  # extract_est_se_ci
  # format_est_se_ci()
  # center_scores()
  # standardize_scores()

  # build_mim_table()
  # build_ao_mim_parameter_table()
  # build_mnsd_table()
  # build_ao_mim_latent_table()
  # build_probit_table()
  # build_prob_at_score_levels_table()

  # hist_item()
  # scatter_p()

# Return vector of mean and standard deviation
mnsd <- function(vals) c(mean(vals, na.rm = TRUE), sd(vals, na.rm = TRUE))

# Rounding that preserves trailing 0s after rounding to d decimal places
rnd <- function(vec, d) format(round(vec, d), nsmall = d)

# Round and format CIs
f_ci <- function(x, d) paste0("[", rnd(x[,1], d), ", ", rnd(x[,2], d), "]")

# Takes in an OpenMx fit object, its summary, and a vector of terms,
# returns a table of estimates, SEs and CI UB and LBs for each parameter
extract_est_se_ci <- function(fitobj, sfit, terms){
  # compute CIs for the coefficients of interest 
  ci <- confint(fitobj)[terms,]
  # swap in the standardized value computed via mxAlgebra (from summary object)
  ci[terms[2],] <- as.numeric(sfit$CI[, c("lbound", "ubound")])
  
  # extract the est and se
  params <- sfit$parameters
  rownames(params) <- params$name
  Est <- params[terms, "Estimate"]
  SE <- params[terms, "Std.Error"]
  
  tb <- cbind(Est, SE, ci)
  tb[terms[2], "Est"] <- sfit$CI[, "estimate"]
  tb[terms[2], "SE"] <- mxSE(beta1, fitobj)
  colnames(tb)[3:4] <- c("LB", "UB")
  return (tb)
}

# Round to the specified number of digits, format CI brackets, and return as 
# a string of alternating Est, SE, CI for each parameter
format_est_se_ci <- function(tb, d) {
  c(rbind(
    rnd(tb[ ,"Est"], d), 
    rnd(tb[ ,"SE"], d), 
    f_ci(tb[, c("LB", "UB")], d)
  ))
}

# Center a vector of scores
center_scores <- function(vec) vec - mean(vec, na.rm = TRUE) 

# Standardize a vector of scores
standardize_scores <- function(vec) {
  (vec - mean(vec, na.rm = TRUE)) / sd(vec, na.rm = TRUE)   
}

# Table functions ####

# Build a table containing information on the configural, partial metric, 
# scalar, and partial strict invariance models tested as part of the traditional
# MIM approach. Round to 'dig' digits.
build_mim_table <- function(dig) {
  tab_fit <- readRDS("rds/tab_fit.rds") # created in `MIM_Step2_detailed.R`
  loadings <- c("i1, i2, i4, i5", "-", "i1", "i2", "i4", "i5",
                "i2, i1","i2, i4", "i2, i5", rep("i2, i1", 12))
  intercepts <- c(rep("i1, i2, i4, i5", 9), rep("i2, i1", 12))
  variances <- c(rep("i1, i2, i4, i5", 10), "-", "i1", "i2", "i4", "i5",
                 "i1, i2", "i1, i4", "i1, i5",  "i1, i2, i4",
                 "i1, i2, i5",  "i1, i2, i4, i5")
  Type <- c("Configural", "Metric", rep("Metric (partial)", 7), "Scalar (partial)", 
            "Strict (partial)", rep("Strict (partial)", 10))
  tab_fit_r <- cbind(rnd(tab_fit[ ,1], dig), tab_fit[, 2:3], rnd(tab_fit[4:8], dig), 
                     tab_fit[, 9], rnd(tab_fit[, 10], dig))
  all_models_table <- cbind(Type, loadings, intercepts, variances,
                            tab_fit_r)
  colnames(all_models_table)[1] <- "Model"
  
  all_models_table[,-1] %>% 
    kbl(booktabs = T, align = "c",#align = "lllcccccccccccc", 
        col.names = c("Loadings" ,"Intercepts", "Uniqueness",  # "$\\lambda$", " $\\nu$", " $\\theta$", 
                      "$\\chi^2$", "df", "npar", "CFI", "TLI", 
                      "RMSEA", "SRMR", "$\\Delta$$\\chi^2$", "$\\Delta$df", "$p$"),
        linesep = "", escape = FALSE) %>%
    pack_rows(index = c("Configural invariance" = 1, "Metric invariance" = 1,
                        "Metric (partial) invariance" = 7, "Scalar (partial) invariance" = 1, 
                        "Strict (partial) invariance" = 1, "Strict (partial) invariance" = 10)) %>%  
    row_spec((c(1, 2, 9, 10, 11)), hline_after = T) %>%
    add_header_above(c("", "Free parameters" = 3, " " = 10)) %>% 
    kable_styling(latex_options = "scale_down") %>%
    kableExtra::footnote(
      general = "The first three columns indicate which parameters were freed on which items in a given model. Cells with '-' indicate that all items were constrained. The final three columns contain the LRT results for comparisons between more and less constrained models (e.g., the LRT for M0 vs. M1 is reported in the second row, while the LRTs between M1 and M2-M5 are reported in rows 3-6). The final partial invariance model is M20. CFI = Comparative Fit Index. TLI = Tucker-Lewis Index. RMSEA = Root Mean Square Error of Approximation. SRMR = Square Root Mean Residual.",
      threeparttable = TRUE, footnote_as_chunk = TRUE) 
  }

# Build a table containing the loading, intercept, and uniqueness estimates 
# from AO and MIM approaches. Round to 'dig' digits.
build_ao_mim_parameter_table <- function(dig) {
  est_partial <- readRDS("rds/est_partial.rds")
  est_align <- readRDS("rds/est_align.rds")
  loadings_tab123 <- rnd(
    cbind(est_align$ELS$lambda, 
          c(est_align$HSLS$lambda[1:2], NA, est_align$HSLS$lambda[3:4]),
          est_partial$ELS$lambda, 
          c(est_partial$HSLS$lambda[1:2], NA, est_partial$HSLS$lambda[3:4])), dig)
  intercepts_tab <- rnd(
    cbind(est_align$ELS$nu, 
          c(est_align$HSLS$nu[1:2], NA, est_align$HSLS$nu[3:4]),
          est_partial$ELS$nu, 
          c(est_partial$HSLS$nu[1:2], NA, 
            est_partial$HSLS$nu[3:4])), dig)
  var_tab123 <- rnd(
    cbind(diag(est_align$ELS$theta), 
          c(diag(est_align$HSLS$theta)[1:2], NA, 
            diag(est_align$HSLS$theta)[3:4]),
          diag(est_partial$ELS$theta), 
          c(diag(est_partial$HSLS$theta)[1:2], NA, 
            diag(est_partial$HSLS$theta)[3:4])), dig)
  
  colnames(intercepts_tab) <- colnames(loadings_tab123) <- colnames(var_tab123) <- NULL 
  binded <- cbind(loadings_tab123[,1:2], intercepts_tab[,1:2], var_tab123[,1:2],
                  loadings_tab123[,3:4], intercepts_tab[,3:4], var_tab123[,3:4])
  colnames(binded) <- NULL
  rownames(binded) <- c("i1", "i2", "i3", "i4", "i5")
  binded2 <- gsub("NA", "-", binded)
  binded_kbl <- 
    kbl(binded2, booktabs = T, align = "c", linesep = "") %>%
    kable_classic() %>%    
    add_header_above(c(" " = 1, "ELS" = 1, "HSLS" = 1, "ELS" = 1, "HSLS" = 1,
                       "ELS" = 1, "HSLS" = 1, "ELS" = 1, "HSLS" = 1, "ELS" = 1, 
                       "HSLS" = 1, "ELS" = 1, "HSLS" = 1)) %>%
    add_header_above(c("  " = 1, "Loadings" = 2, "Intercepts" = 2, "Uniqueness" = 2,
                       "Loadings" = 2, "Intercepts" = 2, "Uniqueness" = 2)) %>%
    add_header_above(c(" " = 1, "AO (Approximate Invariance)" = 6, "MIM (Partial Invariance)" = 6)) %>%
    column_spec(7, border_right = T) %>%
    kableExtra::footnote(general = "An approximate invariance model was determined via alignment optimization (AO) and a partial invariance model was determined via measurement invariance modeling (MIM).",
                         threeparttable = TRUE, footnote_as_chunk = TRUE)
  return(binded_kbl)  
}

# Build a table displaying the means and standard deviations of the standardized
# and unstandardized factor scores from AO, MIM and CMS approaches. Round to 
# 'dig' digits.
build_mnsd_table <- function(capt = NULL, dig) {
  df <- readRDS("rds/score_df_continuous.rds")

   # compute mn, sd of variables in vector v for sample s  
   mnsd_sv <- function(s = NULL, v) {
     if(!is.null(s)) {
       return(unlist(lapply(v, function(x) mnsd(df[df$sample == s, x[1]]))))
     }
     return(unlist(lapply(v, function(x) mnsd(df[, x[1]]))))
   }  
   mnsd_f <- function(v) { # make and format table of mn, sd by sample
      tab <- rnd(rbind(mnsd_sv("ELS", v), mnsd_sv("HSLS", v), mnsd_sv(NULL, v)), dig)
      rownames(tab) <- c("ELS", "HSLS", "overall")
      colnames(tab) <- rep(c("M", "SD"), 3)
      return(tab)
   }
 
  std_v <- c("self_efficacy_ao_c", "self_efficacy_mim_c", "self_efficacy_mn_c")
  raw_v <- c("approx_cont", "partial_cont", "mean_score")
  
  tab <- cbind(mnsd_f(std_v), mnsd_f(raw_v))
  
  tab_kbl <-
    kbl(tab, align = "c", booktabs = TRUE, caption = capt) %>% 
    add_header_above(c("", "AO" = 2, "MIM" = 2, "CMS" = 2, "AO" = 2, "MIM" = 2, "CMS" = 2)) %>%
    add_header_above(c(" " = 1, "Centered" = 6, "Raw" = 6)) %>%
    column_spec(7, border_right = T) %>% kable_classic() %>%
    kableExtra::footnote(general = "Columns labeled AO and MIM refer to math self-efficacy factor scores computed using an approximate invariance model determined via alignment optimization (AO) or a partial invariance model determined via measurement invariance modeling (MIM) respectively. CMS refers to composite mean scores.",
                         threeparttable = TRUE, footnote_as_chunk = TRUE)
  return(tab_kbl)
}

# Build a table containing the latent mean and variances for the AO, MIM
# approaches. Round to 'dig' digits.
build_ao_mim_latent_table <- function(dig) {
  est_partial <- readRDS("rds/est_partial.rds")
  est_align <- readRDS("rds/est_align.rds")
  alpha_psi_tab <- rbind(
    cbind(rnd(est_align$ELS$alpha, dig), rnd(est_align$HSLS$alpha, dig),
          rnd(est_partial$ELS$alpha, dig), rnd(est_partial$HSLS$alpha, dig)),
    cbind(rnd(est_align$ELS$psi, dig), rnd(est_align$HSLS$psi, dig),
          rnd(est_partial$ELS$psi, dig), rnd(est_partial$HSLS$psi, dig)))
  colnames(alpha_psi_tab) <- NULL
  rownames(alpha_psi_tab) <-c("Latent mean", "Latent variance")
  latent_tab_kbl <- 
    kbl(alpha_psi_tab, booktabs = T, align = "c", linesep = "") %>%
    kable_styling(latex_options = "scale_down") %>% 
    kable_classic() %>% 
    add_header_above(c(" " = 1, "ELS" = 1, "HSLS" = 1, "ELS" = 1, "HSLS" = 1)) %>%
    add_header_above(c(" " = 1, "AO (Approximate Invariance)" = 2, "MIM (Partial Invariance)" = 2)) %>% 
    # Equalize column widths
    column_spec(2, width = "2cm") %>%
    column_spec(3, width = "2cm") %>%
    column_spec(4, width = "2cm") %>%
    column_spec(5, width = "2cm") %>%
    kableExtra::footnote(general = "Latent mean and latent variances were constrained to 0 and 1 respectively in ELS. An approximate invariance model was determined via alignment optimization (AO) and a partial invariance model was determined via measurement invariance modeling (MIM).",
                         threeparttable = TRUE, footnote_as_chunk = TRUE)
  return(latent_tab_kbl)  
}

# Build a table displaying the coefficient estimates and fit information for 
# probit regression models predicting dropout from math self-efficacy and sample
# membership
build_probit_table <- function() {
  tab_probit <- readRDS("rds/probit_reg_table_std.rds") # created in subsequent_analyses.R
 # colnames(tab_probit)[5] <- "CMS"
  kbl(tab_probit, align = "c", booktabs = T)  %>%
    row_spec((tab_probit %>% nrow() - 4), hline_after = T) %>%
    pack_rows(index = c("(intercept)" = 3, "self-efficacy"= 3, "ELS"= 3, 
                        "Model information" = 4)) %>%
    add_header_above(c(" ", "FS (PA)" = 2, " " = 1, "FS (2S-PA)" = 2)) %>%
    kable_styling(latex_options = "scale_down") %>% 
    kableExtra::footnote(general = "'self-efficacy' refers to centered Bartlett factor scores (FS) computed using the approximate and partial invariance models determined with the alignment optimization (AO) and measurement invariance modeling (MIM) approaches respectively, or to centered observed composite mean scores (CMS) on the math self-efficacy test items. The coefficients for self-efficacy were standardized, and indicate the change in the probit for a one standard deviation increase in scores. FS were corrected for measurement noninvariance and unreliability following a two-stage path analysis (2S-PA) approach or only for measurement noninvariance (PA). CMS do not account for noninvariance or unreliability. ELS refers to sample membership (ELS: 1, HSLS: 0). Values in this table are reported with three decimal places to highlight the subtle differences across conditions.",
                         threeparttable = TRUE, footnote_as_chunk = TRUE)
}

# Build a table of predicted dropout probabilities by sample at specific math 
# self-efficacy levels
build_prob_at_score_levels_table <- function() {
  tab_std <- readRDS("rds/probit_reg_table_std.rds") # created in subsequent_analyses.R
  
  models <- c("AO", "MIM", "CMS", "AO-2S-PA", "MIM-2S-PA")
  b <- matrix(as.numeric(unlist(tab_std[c(1, 4, 7), models])), ncol = 5) # betas
  

  ao_f <- function(score, ELS) pnorm(b[1, 1] + b[2, 1] * score + b[3, 1] * ELS)
  mim_f <- function(score, ELS) pnorm(b[1, 2] + b[2, 2] * score + b[3, 2] * ELS)
  mn_f <- function(score, ELS) pnorm(b[1, 3] + b[2, 3] * score + b[3, 3] * ELS)
  ao_2spa_f <- function(score, ELS) pnorm(b[1, 4] + b[2, 4] * score + b[3, 4] * ELS)
  mim_2spa_f <- function(score, ELS) pnorm(b[1, 5] + b[2, 5] * score + b[3, 5] * ELS)
  
  vec_p <- function(s, e) {
    df <- cbind(ao_f(s, e), mim_f(s, e), mn_f(s, e),  ao_2spa_f(s, e), mim_2spa_f(s, e))
    colnames(df) <- models
    df * 100
  }
  p_els <- matrix(c(rbind(vec_p(-2, 1), vec_p(-1, 1), vec_p(0, 1), vec_p(1, 1), 
                          vec_p(2, 1))), ncol = 5)
  p_hsls <- matrix(c(rbind(vec_p(-2, 0), vec_p(-1, 0), vec_p(0, 0), vec_p(1, 0), 
                           vec_p(2, 0))), ncol = 5)
  p_diff <- p_hsls - p_els
  p_els_hsls <- matrix(ncol = 0, nrow = nrow(p_els))  # empty matrix
  for (i in 1:ncol(p_els)) {# loop through and interweave columns
    p_els_hsls <- cbind(p_els_hsls, p_els[, i], p_hsls[, i], p_diff[, i])
  }
  p_els_hsls <- cbind(c(-2, -1, 0, 1, 2), p_els_hsls)
  
  rnd(p_els_hsls, 2) %>% 
    kbl(align = "c", booktabs = T, linesep = "", escape = FALSE,
        col.names = c("self-efficacy", rep(c("P(1|E)", "P(1|H)", "$\\Delta$"), 5))) %>% 
    kable_styling(latex_options = "scale_down") %>%
    add_header_above(c(" " = 1, "AO-PA" = 3, "MIM-PA" = 3, "CMS" = 3, "AO-2S-PA" = 3, "MIM-2S-PA" = 3)) %>%
    kableExtra::footnote(general = "Rows represent specific levels of centered math self-efficacy scores. P(1|E) and P(1|H) denote predicted probability of dropout (in percentages) at a given level of math self-efficacy for a student in the ELS and HSLS samples respectively. The difference between these probabilities is given in the third column for each approach. CMS refers to composite mean scores.",
                         threeparttable = TRUE, footnote_as_chunk = TRUE)
}


# Figure functions ####

# Histograms
hist_item <- function(df, color) {
  cols <- c("#F26178", "#FFCC00")
  par(las=2)
  hist(df, main = "", xlab = "Responses", xaxt = "n", col = color, 
       ylim = c(0, 12000), cex.axis = 0.7)#, yaxt = "n")
  axis(1, at = 1:4)
}
# apply(dat[dat$sample == "ELS",m_items], MARGIN = 2,
#       FUN = hist_item, color = cols[1])
# apply(dat[dat$sample == "HSLS", m_items[-3]], MARGIN = 2, 
#       FUN = hist_item, color = cols[2])
# 

# Build a ggplot scatterplot with an optional 45 degree line
# and other specifications
scatter_p <- function(dataset, x, y, xlab, ylab, title, ref_line = TRUE) {
  cols <- c("#F26178", "#FFCC00")
  p <- ggplot(dataset, aes(x, y, col = sample)) +
    xlab(xlab) +  ylab(ylab) + ggtitle(title) + xlim(-2, 2) + ylim(-2,2)+
    theme_classic() + geom_point() + scale_color_manual(values = cols)
  if (ref_line) {
    p <- p + geom_abline(slope = 1, intercept = 0) +
      coord_fixed()
  }
  p + geom_smooth(method = "lm", linewidth = 0.5, se = FALSE)
}

#AO vs. MIM
# scatter_p(df, y = df$self_efficacy_ao_c, x = df$self_efficacy_mim_c,
#           xlab = "math self-efficacy FS (MIM)",
#           ylab = "math self-efficacy FS (AO)",
#           "Correlation of Centered Bartlett Factor Scores (FS) from AO and MIM")
# 


# Notes on syntax:
# linesep = "" disables the automatic line break added after 5 rows for readability
# extra_latex_after = "%" disables the additional line break after hline
# kableExtra::footnote needs to be specified to differentiate from flextable::footnote
# format(round(vec, d), nsmall = d) preserves zeros after the decimal point after
# rounding up
# format = "latex" latex of tables
# column_spec(7, border_right = T) for adding a vertical border
# pack_rows(index = table(binded2[,1])) to group by rows, or add_indent()


# others

# Obtain path estimates, SEs, and CI bounds from umxRAM() output for the 
# specified terms
get_path_ests_umx <- function(terms, umx_out) {
  cols <- c("estimate", "lbound", "ubound")
  est_ci <- summary(umxConfint(umx_out, run = TRUE, parm = "all"))$CI[terms, cols]
  umx_ests <- summary(umx_out)$parameters
  rownames(umx_ests) <- umx_ests$name
  se <- umx_ests[terms, "Std.Error"]
  ests <- cbind("Est" = est_ci[,"estimate"], "SE" = se, 
                "LB" = est_ci[,"lbound"], "UB" = est_ci[,"ubound"])
  rownames(ests) <- c("(Intercept)", "self-efficacy", "ELS")
  ests["(Intercept)", c("Est", "LB", "UB")] <- -ests["(Intercept)", c("Est", "UB", "LB")]
  ests
}

# Obtain path estimates, SEs, and CI bounds from glm() output
get_path_ests_glm <- function(glm_out) {
  cis <- confint(glm_out)
  est_se <- summary(glm_out)$coef[, c("Estimate", "Std. Error")]
  ests <- cbind(est_se, cis)
  colnames(ests) <- c("Est", "SE", "LB", "UB")
  rownames(ests)[3] <- "ELS"
  ests
}

# Get fit indices from umxRAM() output, round to the specified digits
fit_umx <- function(s, d) c(s$numObs, rnd(c(s$AIC, s$BIC, s$Minus2LogLikelihood), d))
