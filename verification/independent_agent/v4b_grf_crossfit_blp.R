# Claim 4 supplement: 5-fold cross-fitted CATEs (forest trained on other folds) then BLP
# Y* AIPW-style: regress (Y - m.hat) on (W - e)*1 and (W - e)*(tau - mean tau), HC robust.
.libPaths(c("~/Rlibs_verifier", .libPaths()))
suppressMessages({library(grf); library(sandwich)})
d <- read.csv("data/raw/dataverse_files/analysis_sample.csv"); d <- d[!is.na(d$hba1c_val_6mo), ]
d$days_a1c_consent <- d$consent_date - d$hba1c_date_0mo
cov <- c("age","hba1c_val_0mo","days_a1c_consent","bmi_val_0mo","sys_val_0mo","dias_val_0mo",
         "op_n_0mo","ed_or_ip_n_0mo","male","hispanic","race_nonhisp_black","rural_site",
         "hba1c_cat_control","drug_insulin_0mo","drug_metformin_0mo","drug_glp1_0mo",
         "drug_sulf_0mo","drug_other_0mo","drug_any_db_0mo","ip_i_0mo","ed_i_0mo")
X <- as.matrix(d[, cov]); for (j in seq_len(ncol(X))) X[is.na(X[, j]), j] <- median(X[, j], na.rm = TRUE)
W <- d$treat_now; e <- ave(W, d$rural_site * 2 + d$hba1c_cat_control); Y <- d$hba1c_val_6mo
out <- c()
for (seed in 1:10) {
  set.seed(seed); f <- sample(rep(1:5, length.out = nrow(X))); tau <- mhat <- rep(NA, nrow(X))
  for (k in 1:5) {
    cf <- causal_forest(X[f != k, ], Y[f != k], W[f != k], W.hat = e[f != k], num.trees = 2000, seed = seed)
    tau[f == k] <- predict(cf, X[f == k, ])$predictions
    rf <- regression_forest(X[f != k, ], Y[f != k], num.trees = 1000, seed = seed)
    mhat[f == k] <- predict(rf, X[f == k, ])$predictions
  }
  A <- W - e; fit <- lm(I(Y - mhat) ~ 0 + A + I(A * (tau - mean(tau))))
  b <- coef(fit)[2]; se <- sqrt(vcovHC(fit, type = "HC3")[2, 2])
  out <- rbind(out, c(seed = seed, slope = b, se = se, p_one = pnorm(b / se, lower.tail = FALSE), sd_tau = sd(tau)))
}
print(round(out, 3)); cat("median slope", median(out[, 2]), " median one-sided p", median(out[, 4]), "\n")
