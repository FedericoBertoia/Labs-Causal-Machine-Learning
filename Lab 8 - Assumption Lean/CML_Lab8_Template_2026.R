#############################################
##  Assumption-Lean Regression
##  King County 2001 Birth Data
#############################################

# ── Packages ───────────────────────────────────────────────────────────────────
library(SuperLearner)
library(lmtest)
library(sandwich)
library(fastDummies)
library(dplyr)
library(ggplot2)
library(tibble)
library(patchwork)

expit <- function(x) 1 / (1 + exp(-x))
logit <- function(p) log(p / (1 - p))

set.seed(123)


#############################################
##  Data Loading & Preparation
#############################################

df <- read.table("KingCounty2001.data")
colnames(df) <- c(
  "gender", "plural", "age", "race", "parity", "married",
  "bwt", "smokeN", "drinkN", "firstep", "welfare",
  "smoker", "drinker", "wpre", "wgain", "education", "gestation"
)

############################################################################

# King County 2001 Birth Data
# =========================================================================
#   
# This data file contains outcomes on 2500 children that were born in
# King County in 2001.  
# 
# The variables are (in column order):
#   
# "gender"        M = male, F = female baby
# "plural"        1 = singleton, 2 = twin, 3 = triplet
# "age"           mother's age in years
# "race"          race categories (for mother)
# "parity"        number of previous live born infants
# "married"       Y = yes, N = no
# "bwt"           birth weight in grams
# "smokeN"        number of cigarettes smoked per day during pregnancy
# "drinkN"        number of alcoholic drinks per week during pregnancy
# "firstep"       1 = participant in program; 0 = did not participate
# "welfare"       1 = participant in public assistance program; 0 = did not 
# "smoker"        Y = yes, N = no, U = unknown
# "drinker"       Y = yes, N = no, U = unknown 
# "wpre"          mother's weight in pounds prior to pregnancy
# "wgain"         mother's weight gain in pounds during pregnancy
# "education"     highest grade completed (add 12 + 1 / year of college)
# "gestation"     weeks from last menses to birth of child

############################################################################

# Binary outcome: low birth weight (< 2500 g)
df$low_bwt <- as.integer(df$bwt < 2500)

# Recode character variables to numeric
df$gender  <- as.integer(df$gender  == "M")   # 1 = male
df$married <- as.integer(df$married == "Y")   # 1 = married

# Dummy-code race (drop reference category "other")
df <- dummy_cols(df, select_columns = "race",
                 remove_first_dummy = FALSE, remove_selected_columns = TRUE)
df$race_other <- NULL

# Drop variables not used as confounders
df <- df %>%
  select(-plural, -bwt, -smokeN, -drinkN, -welfare, -smoker,
         -drinker, -wgain, -gestation, -parity)

N <- nrow(df)

# ── Analysis variables ─────────────────────────────────────────────────────────

# Binary exposure: First Steps program participation 
Y <- df$low_bwt   # outcome
A <- df$firstep   # exposure
L <- df %>% select(-low_bwt, -firstep)   

# Continuous exposure: maternal age 
y <- df$low_bwt   # outcome
a <- df$age       # exposure
l <- df %>% select(-low_bwt, -firstep, -age)   


#############################################
###  Exploratory Data Analysis
#############################################

p1 <- df %>%
  mutate(firstep = factor(firstep, labels = c("Non-participant", "Participant"))) %>%
  ggplot(aes(x = firstep, fill = firstep)) +
  geom_bar(alpha = 0.8, width = 0.4) +
  scale_fill_manual(values = c("steelblue", "tomato")) +
  labs(x = NULL, y = "Count", title = "First Steps participation") +
  theme_minimal() +
  theme(legend.position = "none")

p2 <- df %>%
  mutate(low_bwt = factor(low_bwt, labels = c("Normal weight", "Low birth weight"))) %>%
  ggplot(aes(x = low_bwt, fill = low_bwt)) +
  geom_bar(alpha = 0.8, width = 0.4) +
  scale_fill_manual(values = c("steelblue", "tomato")) +
  labs(x = NULL, y = "Count", title = "Low birth weight outcome") +
  theme_minimal() +
  theme(legend.position = "none")

p1 + p2

ggplot(df, aes(x = age)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30,
                 fill = "steelblue", alpha = 0.7) +
  geom_density(color = "navy", linewidth = 0.8) +
  labs(x = "Maternal age (years)", y = "Density",
       title = "Distribution of maternal age") +
  theme_minimal()

df %>%
  mutate(
    firstep = factor(firstep, labels = c("Non-participant", "Participant")),
    low_bwt = factor(low_bwt, labels = c("Normal", "Low BW"))
  ) %>%
  count(firstep, low_bwt) %>%
  group_by(firstep) %>%
  mutate(prop = n / sum(n)) %>%
  ggplot(aes(x = firstep, y = prop, fill = low_bwt)) +
  geom_col(position = "dodge", alpha = 0.85, width = 0.5) +
  geom_text(aes(label = scales::percent(prop, accuracy = 0.1)),
            position = position_dodge(width = 0.5), vjust = -0.4, size = 3.5) +
  scale_fill_manual(values = c("steelblue", "tomato")) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(x = NULL, y = "Proportion", fill = NULL,
       title = "Low birth weight rate by First Steps participation") +
  theme_minimal()


#############################################
###  SuperLearner Setup & Cross-Fitting Folds
#############################################

SL.library <- c("SL.glm", "SL.ranger", "SL.earth", "SL.gam")

n_folds <- 5
fold_id <- sample(rep(1:n_folds, length.out = N))
folds   <- split(seq_len(N), fold_id)


################################################################################
###  Question 1 — Risk Difference: First Steps → Low Birth Weight
################################################################################

# Which variables do you include as confounders L?
# Are there any variables you must exclude, and why?


# ── Step 1: Cross-fitted nuisance estimates ────────────────────────────────────
#
#   mu_hat(A, L) = E[Y | A, L]   — outcome regression
#   pi_hat(L)   = P(A = 1 | L)  — propensity score

AL <- data.frame(A, L)

mu_hat <- rep(NA_real_, N)   # E[Y | A, L]
pi_hat <- rep(NA_real_, N)   # P(A = 1 | L)


for (k in 1:n_folds) {
  
  test  <- folds[[k]]
  train <- setdiff(seq_len(N), test)
  
  # Outcome model: E[Y | A, L]
  outcome_model <- SuperLearner(
    Y          = ...,
    X          = ...,
    SL.library = SL.library,
    family     = ...,
    method     = "method.NNLS",
    cvControl  = list(V = 5)
  )
  
  # Propensity score model: P(A = 1 | L)
  propensity_model <- SuperLearner(
    Y          = ...,
    X          = ...,
    SL.library = SL.library,
    family     = ...,
    method     = "method.NNLS",
    cvControl  = list(V = 5)
  )
  
  mu_hat[test] <- predict(outcome_model,    newdata = AL[test, ])$pred
  pi_hat[test] <- predict(propensity_model, newdata = L[test, ])$pred
}

# Clip propensity scores to avoid numerical instability
pi_hat <- pmax(0.01, pmin(0.99, pi_hat))
mu_hat <- pmax(0.01, pmin(0.99, mu_hat))


# ── Step 2: Assumption-lean estimator ─────────────────────────────────────────
#
#  Regress residual outcome (Y - mu_hat) on residual exposure (A - pi_hat),
#  no intercept.

# Via lm with sandwich SE
lm_rd1        <- lm(I(...) ~ -1 + I(...))
lm_rd1_result <- coeftest(lm_rd1, vcov = sandwich)

beta_rd1_lm <- lm_rd1_result[1, 1]
se_rd1_lm   <- lm_rd1_result[1, 2]

# Via efficient influence function (EIF)
# Hint: beta = E[(A - pi)(Y - mu)] / E[(A - pi)^2]
beta_rd1_eif <- ...

eif_rd1    <- ...
se_rd1_eif <- sd(eif_rd1) / sqrt(N)

# ── Results ───────────────────────────────────────────────────────────────────
results_q1 <- tibble(
  Method   = c("EIF", "LM (sandwich)"),
  Estimate = c(beta_rd1_eif, beta_rd1_lm),
  SE       = c(se_rd1_eif,   se_rd1_lm),
  CI_Lower = Estimate - 1.96 * SE,
  CI_Upper = Estimate + 1.96 * SE
)

knitr::kable(results_q1, digits = 4,
             caption = "Q1 — Risk difference: First Steps on low birth weight")



################################################################################
###  Question 2 — Overlap-Weighted Baseline Characteristics
################################################################################

# Assumption-lean regression implicitly weights units by w(L) = pi(L)(1 - pi(L)).
# Report weighted summary statistics to characterise the effective study population.

# ── Weighted summary ──────────────────────────────────────────────────────────
# Hint: weighted mean = E[(A - pi_hat)^2 * X] / E[(A - pi_hat)^2]
# Hint: EIF = (A - pi_hat)^2 * (X - w_mean) / E[pi_hat * (1 - pi_hat)]

weighted_summary <- do.call(rbind, lapply(names(L), function(varname) {
  X <- L[[varname]]
  if (!is.numeric(X)) return(NULL)
  
  w_mean <- ...
  eif_w  <- ...
  w_se   <- sd(eif_w) / sqrt(N)
  
  data.frame(
    Variable      = varname,
    Weighted_Mean = round(w_mean, 3),
    SE            = round(w_se,   3),
    CI_Lower      = round(w_mean - 1.96 * w_se, 3),
    CI_Upper      = round(w_mean + 1.96 * w_se, 3)
  )
}))

# ── Classic (unweighted) summary for comparison ───────────────────────────────
classic_summary <- do.call(rbind, lapply(names(L), function(varname) {
  X <- L[[varname]]
  if (!is.numeric(X)) return(NULL)
  
  c_mean <- mean(X)
  c_se   <- sd(X) / sqrt(N)
  
  data.frame(
    Variable     = varname,
    Classic_Mean = round(c_mean, 3),
    SE           = round(c_se,   3),
    CI_Lower     = round(c_mean - 1.96 * c_se, 3),
    CI_Upper     = round(c_mean + 1.96 * c_se, 3)
  )
}))

knitr::kable(weighted_summary, caption = "Q2 — Overlap-weighted baseline characteristics")
knitr::kable(classic_summary,  caption = "Q2 — Unweighted baseline characteristics")


################################################################################
###  Question 3 — Risk Difference: Maternal Age → Low Birth Weight
###
###  Same structure as Q1, but with continuous exposure a = age.
###  pi2_hat(L) = E[a | L] is now a conditional mean, not a propensity score.
################################################################################

# ── Step 1: Cross-fitted nuisance estimates ────────────────────────────────────
#
#   mu2_hat(a, l) = E[y | a, l]   — outcome regression
#   pi2_hat(l)   = E[a | l]       — conditional mean of age (Gaussian)

al <- data.frame(a, l)

mu2_hat <- rep(NA_real_, N)   # E[y | a, l]
pi2_hat <- rep(NA_real_, N)   # E[a | l]
y2_hat <- rep(NA_real_, N)   # E[y | l]

for (k in 1:n_folds) {
  
  test  <- folds[[k]]
  train <- setdiff(seq_len(N), test)
  
  # Outcome model: E[y | a, l]
  outcome_model2 <- SuperLearner(
    Y          = ...,
    X          = ...,
    SL.library = SL.library,
    family     = ...,
    method     = "method.NNLS",
    cvControl  = list(V = 5)
  )
  
  
  outcome_model2_y <- SuperLearner(
    Y          = y[train],
    X          = l[train, ],
    SL.library = SL.library,
    family     = binomial(),
    method     = "method.NNLS",
    cvControl  = list(V = 5)
  )
  
  # Exposure model: E[a | l]
  # Hint: continuous exposure → gaussian family
  exposure_model2 <- SuperLearner(
    Y          = ...,
    X          = ...,
    SL.library = SL.library,
    family     = ...,
    method     = "method.NNLS",
    cvControl  = list(V = 5)
  )
  
  mu2_hat[test] <- predict(outcome_model2,  newdata = al[test, ])$pred
  pi2_hat[test] <- predict(exposure_model2, newdata = l[test, ])$pred
  y2_hat[test] <- predict(outcome_model2_y,  newdata = al[test, ])$pred
}

mu2_hat <- pmax(0.01, pmin(0.99, mu2_hat))

# ── Step 2: Assumption-lean estimator ─────────────────────────────────────────

# Via lm with sandwich SE
lm_rd3        <- lm(I(...) ~ -1 + I(...))
lm_rd3_result <- coeftest(lm_rd3, vcov = sandwich)

beta_rd3_lm <- lm_rd3_result[1, 1]
se_rd3_lm   <- lm_rd3_result[1, 2]

# Via EIF
beta_rd3_eif <- ...

eif_rd3    <- ...
se_rd3_eif <- sd(eif_rd3) / sqrt(N)

# ── Results ───────────────────────────────────────────────────────────────────
results_q3 <- tibble(
  Method   = c("EIF", "LM (sandwich)"),
  Estimate = c(beta_rd3_eif, beta_rd3_lm),
  SE       = c(se_rd3_eif,   se_rd3_lm),
  CI_Lower = Estimate - 1.96 * SE,
  CI_Upper = Estimate + 1.96 * SE
)

knitr::kable(results_q3, digits = 6,
             caption = "Q3 — Risk difference: maternal age on low birth weight")

# Interpret your results as a shift intervention.
# How does the estimand differ from Q1?


################################################################################
###  Question 4 — Risk Ratio Scale
###
###  Model:  log E[Y | A, L] = alpha(L) + beta * A
###
###  The debiased pseudo-outcome is:
###    tilde_Y = log(mu_hat) - q_log_hat + (Y - mu_hat) / mu_hat
###  where q_log_hat(L) = E[log mu_hat(A, L) | L]
###
###  Regress tilde_Y on (A - pi_hat), no intercept → sandwich SE → exponentiate.
################################################################################

# ── First Steps ────────────────────────────────────────────────────────────────

# Step 1: estimate q_log_hat(L) = E[log mu_hat(A, L) | L] via cross-fitting
q_log_hat <- rep(NA_real_, N)

for (k in 1:n_folds) {
  
  test  <- folds[[k]]
  train <- setdiff(seq_len(N), test)
  
  # Re-fit outcome model on training fold
  outcome_model <- SuperLearner(
    Y          = Y[train],
    X          = AL[train, ],
    SL.library = SL.library,
    family     = binomial(),
    method     = "method.NNLS",
    cvControl  = list(V = 5)
  )
  
  log_mu_train <- log(predict(outcome_model, newdata = AL[train, ])$pred)
  
  # Marginalisation model: E[log mu_hat | L]
  qlog_model <- SuperLearner(
    Y          = ...,
    X          = ...,
    SL.library = SL.library,
    family     = gaussian(),
    method     = "method.NNLS",
    cvControl  = list(V = 5)
  )
  
  q_log_hat[test] <- predict(qlog_model, newdata = L[test, ])$pred
}

# Step 2: construct pseudo-outcome and regress
pseudo_log1 <- ...

lm_rr1        <- lm(I(pseudo_log1) ~ -1 + I(A - pi_hat))
lm_rr1_result <- coeftest(lm_rr1, vcov = sandwich)

beta_log1 <- lm_rr1_result[1, 1]
se_log1   <- lm_rr1_result[1, 2]

results_rr1 <- tibble(
  Method   = "LM (sandwich)",
  Estimate = exp(beta_log1),
  SE_log   = se_log1,
  CI_Lower = exp(beta_log1 - 1.96 * se_log1),
  CI_Upper = exp(beta_log1 + 1.96 * se_log1)
)

knitr::kable(results_rr1, digits = 3,
             caption = "Q4 — Risk ratio: First Steps on low birth weight")


# ── Maternal age ──────────────────────────────────────────────────────────────

q2_log_hat <- rep(NA_real_, N)

for (k in 1:n_folds) {
  
  test  <- folds[[k]]
  train <- setdiff(seq_len(N), test)
  
  outcome_model2 <- SuperLearner(
    Y          = y[train],
    X          = al[train, ],
    SL.library = SL.library,
    family     = binomial(),
    method     = "method.NNLS",
    cvControl  = list(V = 5)
  )
  
  log_mu2_train <- log(predict(outcome_model2, newdata = al[train, ])$pred)
  
  qlog_model2 <- SuperLearner(
    Y          = ...,
    X          = ...,
    SL.library = SL.library,
    family     = gaussian(),
    method     = "method.NNLS",
    cvControl  = list(V = 5)
  )
  
  q2_log_hat[test] <- predict(qlog_model2, newdata = l[test, ])$pred
}

pseudo_log2 <- ...

lm_rr2        <- lm(I(pseudo_log2) ~ -1 + I(a - pi2_hat))
lm_rr2_result <- coeftest(lm_rr2, vcov = sandwich)

beta_log2 <- lm_rr2_result[1, 1]
se_log2   <- lm_rr2_result[1, 2]

results_rr2 <- tibble(
  Method   = "LM (sandwich)",
  Estimate = exp(beta_log2),
  SE_log   = se_log2,
  CI_Lower = exp(beta_log2 - 1.96 * se_log2),
  CI_Upper = exp(beta_log2 + 1.96 * se_log2)
)

knitr::kable(results_rr2, digits = 4,
             caption = "Q4 — Risk ratio: maternal age on low birth weight")

# Interpret the risk ratios for both exposures.
# What does RR = 1 correspond to?


################################################################################
###  Question 5 — Odds Ratio Scale
###
###  Model:  logit E[Y | A, L] = alpha(L) + beta * A
###
###  The debiased pseudo-outcome is:
###    tilde_Y = logit(mu_hat) - q_logit_hat + (Y - mu_hat) / (mu_hat * (1 - mu_hat))
###  where q_logit_hat(L) = E[logit mu_hat(A, L) | L]
################################################################################

# ── First Steps ────────────────────────────────────────────────────────────────

# Step 1: estimate q_logit_hat(L) = E[logit mu_hat(A, L) | L]
q_logit_hat <- rep(NA_real_, N)

for (k in 1:n_folds) {
  
  test  <- folds[[k]]
  train <- setdiff(seq_len(N), test)
  
  outcome_model <- SuperLearner(
    Y          = Y[train],
    X          = AL[train, ],
    SL.library = SL.library,
    family     = binomial(),
    method     = "method.NNLS",
    cvControl  = list(V = 5)
  )
  
  logit_mu_train <- logit(predict(outcome_model, newdata = AL[train, ])$pred)
  
  qlogit_model <- SuperLearner(
    Y          = ...,
    X          = ...,
    SL.library = SL.library,
    family     = gaussian(),
    method     = "method.NNLS",
    cvControl  = list(V = 5)
  )
  
  q_logit_hat[test] <- predict(qlogit_model, newdata = L[test, ])$pred
}

# Step 2: construct pseudo-outcome and regress

pseudo_logit1 <- ...

lm_or1        <- lm(I(pseudo_logit1) ~ -1 + I(A - pi_hat))
lm_or1_result <- coeftest(lm_or1, vcov = sandwich)

beta_logit1 <- lm_or1_result[1, 1]
se_logit1   <- lm_or1_result[1, 2]

results_or1 <- tibble(
  Method   = "LM (sandwich)",
  Estimate = exp(beta_logit1),
  SE_log   = se_logit1,
  CI_Lower = exp(beta_logit1 - 1.96 * se_logit1),
  CI_Upper = exp(beta_logit1 + 1.96 * se_logit1)
)

knitr::kable(results_or1, digits = 3,
             caption = "Q5 — Odds ratio: First Steps on low birth weight")


# ── Maternal age ──────────────────────────────────────────────────────────────

q2_logit_hat <- rep(NA_real_, N)

for (k in 1:n_folds) {
  
  test  <- folds[[k]]
  train <- setdiff(seq_len(N), test)
  
  outcome_model2 <- SuperLearner(
    Y          = y[train],
    X          = al[train, ],
    SL.library = SL.library,
    family     = binomial(),
    method     = "method.NNLS",
    cvControl  = list(V = 5)
  )
  
  logit_mu2_train <- logit(predict(outcome_model2, newdata = al[train, ])$pred)
  
  qlogit_model2 <- SuperLearner(
    Y          = ...,
    X          = ...,
    SL.library = SL.library,
    family     = gaussian(),
    method     = "method.NNLS",
    cvControl  = list(V = 5)
  )
  
  q2_logit_hat[test] <- predict(qlogit_model2, newdata = l[test, ])$pred
}

pseudo_logit2 <- ...

lm_or2        <- lm(I(pseudo_logit2) ~ -1 + I(a - pi2_hat))
lm_or2_result <- coeftest(lm_or2, vcov = sandwich)

beta_logit2 <- lm_or2_result[1, 1]
se_logit2   <- lm_or2_result[1, 2]

results_or2 <- tibble(
  Method   = "LM (sandwich)",
  Estimate = exp(beta_logit2),
  SE_log   = se_logit2,
  CI_Lower = exp(beta_logit2 - 1.96 * se_logit2),
  CI_Upper = exp(beta_logit2 + 1.96 * se_logit2)
)

knitr::kable(results_or2, digits = 3,
             caption = "Q5 — Odds ratio: maternal age on low birth weight")

# ── Visualise the OR transformation ───────────────────────────────────────────
# Plot E[Y | A=1, L] as a function of E[Y | A=0, L] implied by beta_logit1.
# What does it mean when the curve nearly coincides with the identity line?

p0_seq <- seq(0.001, 0.999, length.out = 1000)
p1_seq <- expit(logit(p0_seq) + beta_logit1)

plot(p0_seq, p1_seq, type = "l", col = "steelblue", lwd = 2,
     xlab = expression(E(Y ~ "|" ~ A == 0 ~ "," ~ L)),
     ylab = expression(E(Y ~ "|" ~ A == 1 ~ "," ~ L)),
     main = expression("Implied" ~ E(Y ~ "|" ~ A == 1 ~ "," ~ L) ~
                         "= expit(logit(" * hat(mu)[A==0] * ") + " * hat(beta) * ")"))
abline(0, 1, col = "gray50", lty = 2)
legend("topleft", legend = c("OR model", "Identity"),
       col = c("steelblue", "gray50"), lty = c(1, 2), bty = "n")


#############################################
###  Summary of Results
#############################################

summary_all <- tibble(
  Exposure = rep(c("First Steps", "Maternal Age"), each = 3),
  Scale    = rep(c("Risk Difference", "Risk Ratio", "Odds Ratio"), times = 2),
  Estimate = c(beta_rd1_eif,  exp(beta_log1),  exp(beta_logit1),
               beta_rd3_eif,  exp(beta_log2),  exp(beta_logit2)),
  SE       = c(se_rd1_eif,    se_log1,         se_logit1,
               se_rd3_eif,    se_log2,         se_logit2),
  Null     = rep(c(0, 1, 1), 2)
) %>%
  mutate(
    CI_Lower = if_else(Scale == "Risk Difference",
                       Estimate - 1.96 * SE,
                       exp(log(Estimate) - 1.96 * SE)),
    CI_Upper = if_else(Scale == "Risk Difference",
                       Estimate + 1.96 * SE,
                       exp(log(Estimate) + 1.96 * SE))
  )

knitr::kable(
  summary_all %>% select(-Null) %>% mutate(across(where(is.numeric), ~round(.x, 4))),
  caption = "Summary — all assumption-lean regression results"
)

# ── Forest plot ───────────────────────────────────────────────────────────────
summary_all %>%
  mutate(Label = factor(paste(Exposure, "-", Scale),
                        levels = rev(unique(paste(Exposure, "-", Scale))))) %>%
  ggplot(aes(x = Estimate, y = Label, color = Exposure)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0.3, linewidth = 0.8) +
  geom_vline(aes(xintercept = Null), linetype = "dashed", color = "gray50") +
  facet_wrap(~Scale, scales = "free_x", nrow = 1) +
  scale_color_manual(values = c("steelblue", "tomato")) +
  labs(x = "Estimate (95% CI)", y = NULL, color = "Exposure",
       title = "Assumption-Lean Regression — King County 2001",
       subtitle = "Risk difference, risk ratio, and odds ratio scales") +
  theme_minimal() +
  theme(legend.position = "bottom")

