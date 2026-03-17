########################
# 0. SETUP
########################

# Packages
library(twang)
library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(SuperLearner)
library(tmle)
library(sandwich)
library(stdReg)
library(npcausal)

# Reproducibility
set.seed(123)

# Load data
data(lindner)

lindner <- lindner %>%
  mutate(survival_binary = ifelse(sixMonthSurvive == FALSE, 1, 0)) %>%
  select(-sixMonthSurvive, -lifepres)

# Initial inspection
glimpse(lindner)

# Confirm no missing values
stopifnot(all(colSums(is.na(lindner)) == 0))


########################
# 1. EXPLORATORY DATA ANALYSIS
########################

# --- Outcome and treatment marginal distributions ---
p1 <- lindner %>%
  mutate(survival_binary = factor(survival_binary, labels = c("Survived", "Died"))) %>%
  ggplot(aes(x = survival_binary, fill = survival_binary)) +
  geom_bar(alpha = 0.8, width = 0.4) +
  scale_fill_manual(values = c("steelblue", "tomato")) +
  theme_minimal() +
  labs(x = NULL, y = "Count", fill = "Outcome", title = "Outcome")

p2 <- lindner %>%
  mutate(abcix = factor(abcix, labels = c("No abciximab", "Abciximab"))) %>%
  ggplot(aes(x = abcix, fill = abcix)) +
  geom_bar(alpha = 0.8, width = 0.4) +
  scale_fill_manual(values = c("steelblue", "tomato")) +
  theme_minimal() +
  labs(x = NULL, y = "Count", fill = "Treatment", title = "Treatment")

options(device = "windows") #options(device = "quartz") on mac
p1 + p2


# --- Survival by treatment group ---
lindner %>%
  mutate(
    abcix           = factor(abcix, labels = c("No abciximab", "Abciximab")),
    survival_binary = factor(survival_binary, labels = c("Survived", "Died"))
  ) %>%
  count(abcix, survival_binary) %>%
  group_by(abcix) %>%
  mutate(prop = n / sum(n)) %>%
  ggplot(aes(x = abcix, y = prop, fill = survival_binary)) +
  geom_col(position = "dodge", alpha = 0.8, width = 0.5) +
  geom_text(aes(label = scales::percent(prop, accuracy = 0.1)),
            position = position_dodge(width = 0.5), vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = c("steelblue", "tomato")) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.05)) +
  theme_minimal() +
  labs(x = NULL, y = "Proportion", fill = "Outcome",
       title = "Survival at 6 months by treatment group")

# --- Covariate distributions by treatment ---
cat_vars <- c("stent", "female", "diabetic", "acutemi", "ves1proc")

lindner %>%
  select(-survival_binary) %>%
  mutate(abcix = factor(abcix, labels = c("No abciximab", "Abciximab"))) %>%
  pivot_longer(cols = -abcix, names_to = "variable", values_to = "value") %>%
  mutate(is_cat = variable %in% cat_vars) %>%
  ggplot(aes(x = value, fill = abcix, color = abcix)) +
  geom_bar(
    data  = ~ filter(.x, is_cat),
    fill  = "steelblue", color = "steelblue", alpha = 0.6,
    width = 0.4
  ) +
  geom_density(
    data  = ~ filter(.x, !is_cat),
    fill  = "steelblue", color = "steelblue", alpha = 0.3, linewidth = 0.8
  ) +
  facet_wrap(~variable, scales = "free") +
  theme_minimal() +
  labs(x = NULL, y = "Count / Density")


########################
# 2. UNADJUSTED ESTIMATORS
########################

# Model formula (unadjusted)
unadjusted <- survival_binary ~ abcix

# --- Risk Difference (linear probability model) ---
lpm_fit <- glm(unadjusted, data = lindner, family = gaussian())
summary(lpm_fit)
confint(lpm_fit)

rd    <- coef(lpm_fit)["abcix"]
ci_rd <- confint(lpm_fit)["abcix", ]

cat("Risk Difference:", round(rd, 4), "\n")
cat("95% CI: [", round(ci_rd[1], 4), ",", round(ci_rd[2], 4), "]\n")

# --- Odds Ratio (logistic regression) ---
logit_fit <- glm(unadjusted, data = lindner, family = binomial())
summary(logit_fit)

or    <- exp(coef(logit_fit)["abcix"])
ci_or <- exp(confint(logit_fit)["abcix", ])

cat("Odds Ratio:", round(or, 4), "\n")
cat("95% CI: [", round(ci_or[1], 4), ",", round(ci_or[2], 4), "]\n")


########################
# 3. SETUP FOR CAUSAL ESTIMATORS
########################

# Outcome, treatment, covariates
Y <- lindner$survival_binary
A <- lindner$abcix
X <- lindner %>% select(-survival_binary, -abcix)

N <- nrow(lindner)

# Candidate learners
SL.library <- c("SL.mean", "SL.glm", "SL.ranger", "SL.gam", "SL.earth")

# Cross-fitting folds
n_folds <- 5
fold_id <- sample(rep(1:n_folds, length.out = N))
folds   <- split(seq_len(N), fold_id)


########################
# 4. G-COMPUTATION
########################

# ── 4.1 Parametric outcome model ──────────────────────────────────────────────

# Model formula
formula_gcomp <- as.formula(paste("Y ~ A +", paste(colnames(X), collapse = " + ")))

# Fit parametric outcome model
glm_Q <- glm(formula_gcomp, family = binomial(), data = data.frame(Y, A, X))

# Predict potential outcomes
dat0      <- data.frame(Y, A = 0, X)
dat1      <- data.frame(Y, A = 1, X)
Q0_par    <- predict(glm_Q, newdata = dat0, type = "response")
Q1_par    <- predict(glm_Q, newdata = dat1, type = "response")

ate_gcomp_par <- mean(Q1_par - Q0_par)
cat("G-Computation ATE (parametric, by hand):", round(ate_gcomp_par, 4), "\n")

# Standard error via stdReg (sandwich estimator)
fit_std <- stdGlm(fit = glm_Q, data = data.frame(Y, A, X), X = "A")
print(summary(fit_std))

ate_gcomp_par_stdreg <- fit_std$est[2] - fit_std$est[1]
ate_gcomp_par_se     <- sqrt(fit_std$vcov[1] + fit_std$vcov[4] - 2 * fit_std$vcov[3])
ci_gcomp_par         <- ate_gcomp_par + c(-1, 1) * 1.96 * ate_gcomp_par_se

cat("G-Computation ATE (parametric, stdReg):", round(ate_gcomp_par_stdreg, 4), "\n")
cat("SE:", round(ate_gcomp_par_se, 4), "\n")
cat("95% CI: [", round(ci_gcomp_par[1], 4), ",", round(ci_gcomp_par[2], 4), "]\n")

# ── 4.2 Non-parametric outcome model (Super Learner) ─────────────────────────

# Without cross-fitting
sl_Q <- SuperLearner(
  Y          = Y,
  X          = data.frame(A = A, X),
  SL.library = SL.library,
  family     = binomial(),
  method     = "method.NNLS",
  cvControl  = list(V = 5)
)

Q1 <- predict(sl_Q, newdata = data.frame(A = 1, X))$pred
Q0 <- predict(sl_Q, newdata = data.frame(A = 0, X))$pred

ate_gcomp <- mean(Q1 - Q0)
cat("G-Computation ATE (no cross-fitting):", round(ate_gcomp, 4), "\n")

# With cross-fitting
Q1_cf <- numeric(N)
Q0_cf <- numeric(N)

for (k in 1:n_folds) {
  
  test  <- folds[[k]]
  train <- setdiff(seq_len(N), test)
  
  sl_Q_cf <- SuperLearner(
    Y          = Y[train],
    X          = data.frame(A = A[train], X[train, ]),
    SL.library = SL.library,
    family     = binomial(),
    method     = "method.NNLS",
    cvControl  = list(V = 5)
  )
  
  Q1_cf[test] <- predict(sl_Q_cf, newdata = data.frame(A = 1, X[test, ]))$pred
  Q0_cf[test] <- predict(sl_Q_cf, newdata = data.frame(A = 0, X[test, ]))$pred
}

ate_gcomp_cf <- mean(Q1_cf - Q0_cf)
cat("G-Computation ATE (cross-fitting):", round(ate_gcomp_cf, 4), "\n")


########################
# 5. IPW
########################

# ── 5.1 Parametric propensity score model ─────────────────────────────────────

# Fit parametric propensity score model
formula_pi  <- as.formula(paste("A ~", paste(colnames(X), collapse = " + ")))
glm_pi      <- glm(formula_pi, family = binomial(), data = data.frame(A, X))
pi_hat_par  <- glm_pi$fitted.values

# Positivity check
summary(pi_hat_par)

data.frame(pi_hat = pi_hat_par, A = factor(A, labels = c("No abciximab", "Abciximab"))) %>%
  ggplot(aes(x = pi_hat, fill = A, color = A)) +
  geom_density(alpha = 0.3, linewidth = 0.8) +
  scale_fill_manual(values  = c("steelblue", "tomato")) +
  scale_color_manual(values = c("steelblue", "tomato")) +
  theme_minimal() +
  labs(x     = expression(hat(pi)(X)),
       y     = "Density",
       fill  = "Treatment",
       color = "Treatment",
       title = "Propensity scores by treatment group (parametric)")

data.frame(
  weight = c(1 / pi_hat_par[A == 1], 1 / (1 - pi_hat_par[A == 0])),
  A      = factor(c(rep(1, sum(A == 1)), rep(0, sum(A == 0))),
                  labels = c("No abciximab", "Abciximab"))
) %>%
  ggplot(aes(x = A, y = weight, fill = A)) +
  geom_boxplot(alpha = 0.6, outlier.shape = 21) +
  scale_fill_manual(values = c("steelblue", "tomato")) +
  theme_minimal() +
  labs(x = NULL, y = "Inverse weight", fill = "Treatment",
       title = "Distribution of inverse weights by treatment group (parametric)")

# Weighted linear regression (Hajek / self-normalized IPW)
weights_par    <- A / pi_hat_par + (1 - A) / (1 - pi_hat_par)
fit_ipw_par    <- lm(Y ~ A, weights = weights_par)
ate_ipw_par    <- coef(fit_ipw_par)["A"]
ate_ipw_par_se <- sqrt(vcovHC(fit_ipw_par, type = "HC")["A", "A"])
ci_ipw_par     <- ate_ipw_par + c(-1, 1) * 1.96 * ate_ipw_par_se

cat("IPW ATE (parametric, weighted regression):", round(ate_ipw_par, 4), "\n")
cat("SE:", round(ate_ipw_par_se, 4), "\n")
cat("95% CI: [", round(ci_ipw_par[1], 4), ",", round(ci_ipw_par[2], 4), "]\n")

# ── 5.2 Non-parametric propensity score model (Super Learner) ─────────────────

# Without cross-fitting
sl_pi <- SuperLearner(
  Y          = A,
  X          = X,
  SL.library = SL.library,
  family     = binomial(),
  method     = "method.NNLS",
  cvControl  = list(V = 5)
)

pi_hat <- predict(sl_pi, newdata = X)$pred

# Positivity check
summary(pi_hat)

data.frame(pi_hat = pi_hat, A = factor(A, labels = c("No abciximab", "Abciximab"))) %>%
  ggplot(aes(x = pi_hat, fill = A, color = A)) +
  geom_density(alpha = 0.3, linewidth = 0.8) +
  scale_fill_manual(values  = c("steelblue", "tomato")) +
  scale_color_manual(values = c("steelblue", "tomato")) +
  theme_minimal() +
  labs(x     = expression(hat(pi)(X)),
       y     = "Density",
       fill  = "Treatment",
       color = "Treatment",
       title = "Propensity scores by treatment group (no cross-fitting)")

data.frame(
  weight = c(1 / pi_hat[A == 1], 1 / (1 - pi_hat[A == 0])),
  A      = factor(c(rep(1, sum(A == 1)), rep(0, sum(A == 0))),
                  labels = c("No abciximab", "Abciximab"))
) %>%
  ggplot(aes(x = A, y = weight, fill = A)) +
  geom_boxplot(alpha = 0.6, outlier.shape = 21) +
  scale_fill_manual(values = c("steelblue", "tomato")) +
  theme_minimal() +
  labs(x = NULL, y = "Inverse weight", fill = "Treatment",
       title = "Distribution of inverse weights by treatment group (no cross-fitting)")

# Horvitz-Thompson estimator
ipw_weights <- A / pi_hat - (1 - A) / (1 - pi_hat)
ate_ipw     <- mean(ipw_weights * Y)
cat("IPW ATE - Horvitz-Thompson (no cross-fitting):", round(ate_ipw, 4), "\n")

# Weighted linear regression
weights_sl  <- A / pi_hat + (1 - A) / (1 - pi_hat)
fit_ipw_sl  <- lm(Y ~ A, weights = weights_sl)

ate_ipw_sl    <- coef(fit_ipw_sl)["A"]
ate_ipw_sl_se <- sqrt(vcovHC(fit_ipw_sl, type = "HC")["A", "A"])
ci_ipw_sl     <- ate_ipw_sl + c(-1, 1) * 1.96 * ate_ipw_sl_se

cat("IPW ATE - weighted regression (no cross-fitting):", round(ate_ipw_sl, 4), "\n")
cat("SE:", round(ate_ipw_sl_se, 4), "\n")
cat("95% CI: [", round(ci_ipw_sl[1], 4), ",", round(ci_ipw_sl[2], 4), "]\n")

# With cross-fitting
pi_hat_cf <- numeric(N)

for (k in 1:n_folds) {
  
  test  <- folds[[k]]
  train <- setdiff(seq_len(N), test)
  
  sl_pi_cf <- SuperLearner(
    Y          = A[train],
    X          = X[train, ],
    SL.library = SL.library,
    family     = binomial(),
    method     = "method.NNLS",
    cvControl  = list(V = 5)
  )
  
  pi_hat_cf[test] <- predict(sl_pi_cf, newdata = X[test, ])$pred
}

# Positivity check
summary(pi_hat_cf)

data.frame(pi_hat_cf = pi_hat_cf, A = factor(A, labels = c("No abciximab", "Abciximab"))) %>%
  ggplot(aes(x = pi_hat_cf, fill = A, color = A)) +
  geom_density(alpha = 0.3, linewidth = 0.8) +
  scale_fill_manual(values  = c("steelblue", "tomato")) +
  scale_color_manual(values = c("steelblue", "tomato")) +
  theme_minimal() +
  labs(x     = expression(hat(pi)(X)),
       y     = "Density",
       fill  = "Treatment",
       color = "Treatment",
       title = "Propensity scores by treatment group (cross-fitting)")

data.frame(
  weight = c(1 / pi_hat_cf[A == 1], 1 / (1 - pi_hat_cf[A == 0])),
  A      = factor(c(rep(1, sum(A == 1)), rep(0, sum(A == 0))),
                  labels = c("No abciximab", "Abciximab"))
) %>%
  ggplot(aes(x = A, y = weight, fill = A)) +
  geom_boxplot(alpha = 0.6, outlier.shape = 21) +
  scale_fill_manual(values = c("steelblue", "tomato")) +
  theme_minimal() +
  labs(x = NULL, y = "Inverse weight", fill = "Treatment",
       title = "Distribution of inverse weights by treatment group (cross-fitting)")

# Horvitz-Thompson estimator
ipw_weights_cf <- A / pi_hat_cf - (1 - A) / (1 - pi_hat_cf)
ate_ipw_cf     <- mean(ipw_weights_cf * Y)
cat("IPW ATE - Horvitz-Thompson (cross-fitting):", round(ate_ipw_cf, 4), "\n")

# Weighted linear regression
weights_sl_cf  <- A / pi_hat_cf + (1 - A) / (1 - pi_hat_cf)
fit_ipw_sl_cf  <- lm(Y ~ A, weights = weights_sl_cf)

ate_ipw_cf_sl    <- coef(fit_ipw_sl_cf)["A"]
ate_ipw_cf_sl_se <- sqrt(vcovHC(fit_ipw_sl_cf, type = "HC")["A", "A"])
ci_ipw_cf_sl     <- ate_ipw_cf_sl + c(-1, 1) * 1.96 * ate_ipw_cf_sl_se

cat("IPW ATE - weighted regression (cross-fitting):", round(ate_ipw_cf_sl, 4), "\n")
cat("SE:", round(ate_ipw_cf_sl_se, 4), "\n")
cat("95% CI: [", round(ci_ipw_cf_sl[1], 4), ",", round(ci_ipw_cf_sl[2], 4), "]\n")


########################
# 6. AIPW
########################

# ── 6.1 npcausal AIPW ────────────────────────────────────────────────────────
# install.packages("devtools")
# library(devtools)
# install_github("ehkennedy/npcausal")

ate_npcausal <- ate(
  y       = Y,
  a       = A,
  x       = X,
  nsplits = 5,
  sl.lib  = SL.library
)


# ── 6.2 Manual AIPW ──────────────────────────────────────────────────────────

# Without cross-fitting
aug_1    <-  A      * (Y - Q1) / pi_hat
aug_0    <- (1 - A) * (Y - Q0) / (1 - pi_hat)
psi_i    <- (Q1 - Q0) + aug_1 - aug_0
ate_aipw <- mean(psi_i)
se_aipw  <- sd(psi_i) / sqrt(N)
ci_aipw  <- ate_aipw + c(-1, 1) * 1.96 * se_aipw

cat("AIPW ATE (no cross-fitting):", round(ate_aipw, 4), "\n")
cat("95% CI: [", round(ci_aipw[1], 4), ",", round(ci_aipw[2], 4), "]\n")

# With cross-fitting
aug_1_cf    <-  A      * (Y - Q1_cf) / pi_hat_cf
aug_0_cf    <- (1 - A) * (Y - Q0_cf) / (1 - pi_hat_cf)
psi_i_cf    <- (Q1_cf - Q0_cf) + aug_1_cf - aug_0_cf
ate_aipw_cf <- mean(psi_i_cf)
se_aipw_cf  <- sd(psi_i_cf) / sqrt(N)
ci_aipw_cf  <- ate_aipw_cf + c(-1, 1) * 1.96 * se_aipw_cf

cat("AIPW ATE (cross-fitting):", round(ate_aipw_cf, 4), "\n")
cat("95% CI: [", round(ci_aipw_cf[1], 4), ",", round(ci_aipw_cf[2], 4), "]\n")


########################
# 7. TMLE
########################

# ── 7.1 tmle package ─────────────────────────────────────────────────────────

# Cross fitting estimating again the nuisances
tmle_fit <- tmle(
  Y            = Y,
  A            = A,
  W            = as.data.frame(X),
  Q.SL.library = SL.library,
  g.SL.library = SL.library,
  family       = "binomial",
  V.Q          = 5,
  V.g          = 5
)

ate_tmle <- tmle_fit$estimates$ATE$psi
ci_tmle  <- tmle_fit$estimates$ATE$CI

cat("TMLE ATE:", round(ate_tmle, 4), "\n")
cat("95% CI: [", round(ci_tmle[1], 4), ",", round(ci_tmle[2], 4), "]\n")

# With cross-fitting (supply cross-fitted nuisance estimates directly)
tmle_fit_cf <- tmle(
  Y            = Y,
  A            = A,
  W            = as.data.frame(X),
  Q.SL.library = SL.library,
  g.SL.library = SL.library,
  family       = "binomial",
  Q            = cbind(Q0_cf, Q1_cf),
  g1W          = pi_hat_cf
)

ate_tmle_cf <- tmle_fit_cf$estimates$ATE$psi
ci_tmle_cf  <- tmle_fit_cf$estimates$ATE$CI

cat("TMLE ATE (cross-fitting):", round(ate_tmle_cf, 4), "\n")
cat("95% CI: [", round(ci_tmle_cf[1], 4), ",", round(ci_tmle_cf[2], 4), "]\n")

# ── 7.2 Manual TMLE ──────────────────────────────────────────────────────────

# Clever covariates
H1 <- A / pi_hat_cf
H0 <- (1 - A) / (1 - pi_hat_cf)

# Estimate fluctuation parameters (one per arm)
delta1 <- coef(glm(Y ~ -1 + offset(qlogis(Q1_cf)) + H1, family = binomial()))
delta0 <- coef(glm(Y ~ -1 + offset(qlogis(Q0_cf)) + H0, family = binomial()))

# Update initial predictions
Q1_cf_updated <- plogis(qlogis(Q1_cf) + delta1 / pi_hat_cf)
Q0_cf_updated <- plogis(qlogis(Q0_cf) + delta0 / (1 - pi_hat_cf))

# ATE estimate
ate_tmle_cf_manual <- mean(Q1_cf_updated - Q0_cf_updated)

# Efficient influence function and standard error
eif <- (Q1_cf_updated - Q0_cf_updated - ate_tmle_cf_manual) +
  A / pi_hat_cf * (Y - Q1_cf_updated) -
  (1 - A) / (1 - pi_hat_cf) * (Y - Q0_cf_updated)

se_tmle_cf_manual <- sd(eif) / sqrt(N)
ci_tmle_cf_manual <- ate_tmle_cf_manual + c(-1, 1) * 1.96 * se_tmle_cf_manual

cat("TMLE ATE (cross-fitting, manual):", round(ate_tmle_cf_manual, 4), "\n")
cat("SE:", round(se_tmle_cf_manual, 4), "\n")
cat("95% CI: [", round(ci_tmle_cf_manual[1], 4), ",", round(ci_tmle_cf_manual[2], 4), "]\n")


########################
# 8. COMPARISON OF ESTIMATORS
########################

results <- tibble(
  Estimator     = c("G-Computation", "G-Computation",
                    "IPW",           "IPW",
                    "AIPW",          "AIPW",
                    "TMLE",          "TMLE"),
  Cross_fitting = rep(c("No", "Yes"), 4),
  ATE           = c(ate_gcomp,    ate_gcomp_cf,
                    ate_ipw,      ate_ipw_cf,
                    ate_aipw,     ate_aipw_cf,
                    ate_tmle,     ate_tmle_cf),
  CI_lower      = c(NA, NA, NA, NA,
                    ci_aipw[1],    ci_aipw_cf[1],
                    ci_tmle[1],    ci_tmle_cf[1]),
  CI_upper      = c(NA, NA, NA, NA,
                    ci_aipw[2],    ci_aipw_cf[2],
                    ci_tmle[2],    ci_tmle_cf[2])
)

print(results)

# Forest-style comparison plot
results %>%
  mutate(Estimator = factor(Estimator, levels = rev(c("G-Computation", "IPW", "AIPW", "TMLE")))) %>%
  ggplot(aes(x = ATE, y = Estimator, color = Cross_fitting, shape = Cross_fitting)) +
  geom_point(size = 3, position = position_dodge(width = 0.4)) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper),
                 height = 0.2, linewidth = 0.8,
                 position = position_dodge(width = 0.4),
                 na.rm = TRUE) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("steelblue", "tomato")) +
  theme_minimal() +
  labs(x = "Estimated ATE", y = NULL, color = "Cross-fitting", shape = "Cross-fitting",
       title = "ATE estimates across estimators and cross-fitting strategies")


########################
# BONUS: AIPW AND TMLE FOR ATT
########################

# ── AIPW for ATT (npcausal) ───────────────────────────────────────────────────
att_npcausal <- att(
  y       = Y,
  a       = A,
  x       = X,
  nsplits = 5,
  sl.lib  = SL.library
)


# ── Manual AIPW for ATT ───────────────────────────────────────────────────────
att_cf_manual <- mean(
  (A / mean(A)) * (Y - Q0_cf) -
    ((1 - A) * pi_hat_cf / (1 - pi_hat_cf)) * (Y - Q0_cf) / mean(A)
)

# Influence function and standard error
att_ic <- (A / mean(A)) * (Y - Q0_cf - att_cf_manual) -
  ((1 - A) * pi_hat_cf / (1 - pi_hat_cf)) * (Y - Q0_cf) / mean(A)

se_att_cf_manual <- sd(att_ic) / sqrt(N)

cat("AIPW ATT (cross-fitting, manual):", round(att_cf_manual, 4), "\n")
cat("SE:", round(se_att_cf_manual, 4), "\n")

# ── TMLE for ATT (via tmle package) ──────────────────────────────────────────
att_tmle_cf    <- tmle_fit_cf$estimates$ATT$psi
ci_att_tmle_cf <- tmle_fit_cf$estimates$ATT$CI

cat("TMLE ATT (cross-fitting):", round(att_tmle_cf, 4), "\n")
cat("95% CI: [", round(ci_att_tmle_cf[1], 4), ",", round(ci_att_tmle_cf[2], 4), "]\n")

# ── Manual TMLE for ATT ───────────────────────────────────────────────────────

# Fluctuation model for ATT (linear model on control arm only)
fluctmod_att <- lm(
  Y[A == 0] ~ -1 + offset(Q0_cf[A == 0]) + I(pi_hat_cf[A == 0] / (1 - pi_hat_cf[A == 0]))
)

# Update initial predictions
Q0_cf_updated_att <- Q0_cf + coef(fluctmod_att) * pi_hat_cf / (1 - pi_hat_cf)

# ATT estimate
att_tmle_cf_manual <- mean((A / mean(A)) * (Y - Q0_cf_updated_att))

# Influence function and standard error
att_ic_tmle <- (A / mean(A)) * (Y - Q0_cf_updated_att - att_tmle_cf_manual) -
  ((1 - A) * pi_hat_cf / (1 - pi_hat_cf)) * (Y - Q0_cf_updated_att) / mean(A)

se_att_tmle_cf_manual <- sd(att_ic_tmle) / sqrt(N)

cat("TMLE ATT (cross-fitting, manual):", round(att_tmle_cf_manual, 4), "\n")
cat("SE:", round(se_att_tmle_cf_manual, 4), "\n")