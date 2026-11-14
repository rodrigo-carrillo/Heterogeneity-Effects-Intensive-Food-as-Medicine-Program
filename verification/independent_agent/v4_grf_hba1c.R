# Claim 4: causal forest on HbA1c (6 mo), independent implementation with R grf.
# Honest OOB predictions; known stratum propensity supplied as W.hat; clusters none.
.libPaths(c("~/Rlibs_verifier", .libPaths()))
suppressMessages(library(grf))
d <- read.csv("data/raw/dataverse_files/analysis_sample.csv")
d <- d[!is.na(d$hba1c_val_6mo), ]
d$days_a1c_consent <- d$consent_date - d$hba1c_date_0mo
cov <- c("age","hba1c_val_0mo","days_a1c_consent","bmi_val_0mo","sys_val_0mo","dias_val_0mo",
         "op_n_0mo","ed_or_ip_n_0mo","male","hispanic","race_nonhisp_black","rural_site",
         "hba1c_cat_control","drug_insulin_0mo","drug_metformin_0mo","drug_glp1_0mo",
         "drug_sulf_0mo","drug_other_0mo","drug_any_db_0mo","ip_i_0mo","ed_i_0mo")
X <- as.matrix(d[, cov]); for (j in seq_len(ncol(X))) X[is.na(X[, j]), j] <- median(X[, j], na.rm = TRUE)
W <- d$treat_now
strat <- d$rural_site * 2 + d$hba1c_cat_control
What <- ave(W, strat)  # known stratum propensity
cat("N =", nrow(d), " treated =", sum(W), " missing filled:", sum(is.na(d[, cov])), "\n")

run <- function(Y, label, seed) {
  set.seed(seed)
  cf <- causal_forest(X, Y, W, W.hat = What, num.trees = 4000, seed = seed)
  ate <- average_treatment_effect(cf)
  cal <- test_calibration(cf)
  tau <- predict(cf)$predictions  # OOB
  # RATE with cross-fitting: split sample, train priority forest on one half, evaluate on other.
  idx <- sample(nrow(X)); h1 <- idx[1:floor(nrow(X)/2)]; h2 <- setdiff(seq_len(nrow(X)), h1)
  rp <- c()
  for (s in list(list(a = h1, b = h2), list(a = h2, b = h1))) {
    cfa <- causal_forest(X[s$a, ], Y[s$a], W[s$a], W.hat = What[s$a], num.trees = 2000, seed = seed)
    cfb <- causal_forest(X[s$b, ], Y[s$b], W[s$b], W.hat = What[s$b], num.trees = 2000, seed = seed)
    pr <- predict(cfa, X[s$b, ])$predictions
    r <- rank_average_treatment_effect(cfb, pr, target = "AUTOC")
    rp <- rbind(rp, c(r$estimate, r$std.err))
  }
  # Combine the two halves (average estimate; SE assuming independence)
  est <- mean(rp[, 1]); se <- sqrt(sum(rp[, 2]^2)) / 2
  # Also OOB-priority RATE on full forest (common, though mildly optimistic)
  r_oob <- rank_average_treatment_effect(cf, tau, target = "AUTOC")
  cat(sprintf("\n== %s (seed %d)\nATE (AIPW) %.3f (se %.3f)\n", label, seed, ate[1], ate[2]))
  print(cal)
  b <- cal[2, ]
  cat(sprintf("diff.forest slope %.3f, one-sided p = %.3f\n", b[1], b[4]))
  cat(sprintf("RATE AUTOC cross-fit: est %.3f se %.3f one-sided p = %.3f\n", est, se, pnorm(est / se, lower.tail = FALSE)))
  cat(sprintf("RATE AUTOC OOB-priority: est %.3f se %.3f one-sided p = %.3f\n",
              r_oob$estimate, r_oob$std.err, pnorm(r_oob$estimate / r_oob$std.err, lower.tail = FALSE)))
  cat(sprintf("OOB tau: sd %.3f, range [%.2f, %.2f]\n", sd(tau), min(tau), max(tau)))
  c(slope = b[1], p_cal = b[4], p_rate = pnorm(est / se, lower.tail = FALSE))
}
res <- rbind(
  run(d$hba1c_val_6mo, "HbA1c level at 6 mo", 1),
  run(d$hba1c_val_6mo, "HbA1c level at 6 mo", 2),
  run(d$hba1c_val_6mo, "HbA1c level at 6 mo", 3),
  run(d$hba1c_val_6mo - d$hba1c_val_0mo, "HbA1c change", 1),
  run(d$hba1c_val_6mo - d$hba1c_val_0mo, "HbA1c change", 2))
print(round(res, 3))
