#############################################
###  Assumption-Lean Regression
###  King County 2001 Birth Data
#############################################

# ── Packages ──────────────────────────────────────────────────────────────────
library(SuperLearner)  # install.packages("SuperLearner")
library(lmtest)        # install.packages("lmtest")
library(sandwich)      # install.packages("sandwich")
library(fastDummies)
library(dplyr)
library(ggplot2)
library(tibble)
library(patchwork)

# ── Reproducibility ───────────────────────────────────────────────────────────
set.seed(123)

# ── Load data ─────────────────────────────────────────────────────────────────
# Update the path below to point to your local copy of KingCounty2001.data
file_path <- "KingCounty2001.data"

df <- read.table(file_path)
colnames(df) <- c("gender", "plural", "age", "race", "parity", "married",
                  "bwt", "smokeN", "drinkN", "firstep", "welfare",
                  "smoker", "drinker", "wpre", "wgain", "education", "gestation")

# ── Initial inspection ────────────────────────────────────────────────────────
glimpse(df)


#############################################
###  Data Preparation
#############################################

N <- nrow(df)

# Binary outcome: low birth weight (< 2500 g)
df$low_bwt <- as.integer(df$bwt < 2500)

# Recode character variables to numeric
df$gender  <- as.integer(df$gender  == "M")   # 1 = male
df$married <- as.integer(df$married == "Y")   # 1 = married
df <- dummy_cols(df,
                 select_columns          = "race",
                 remove_first_dummy      = FALSE,
                 remove_selected_columns = TRUE)

df$race_other <- NULL

# Drop unused columns from df permanently
df <- df %>% select(-plural, -bwt, -smokeN, -drinkN, -welfare, -smoker, -drinker, -wgain, -gestation)

# Q1 & Q2: exposure = First Steps (binary), outcome = low_bwt
Y <- df$low_bwt
A <- df$firstep
L <- df %>% select(-firstep, -low_bwt)

# Q3: exposure = maternal age (continuous), outcome = low_bwt
y <- df$low_bwt
a <- df$age
l <- df %>% select(-firstep, -low_bwt, -age) 


#############################################
###  Exploratory Data Analysis
#############################################

p1 <- df %>%
  mutate(firstep = factor(firstep, labels = c("Non-participant", "Participant"))) %>%
  ggplot(aes(x = firstep, fill = firstep)) +
  geom_bar(alpha = 0.8, width = 0.4) +
  scale_fill_manual(values = c("steelblue", "tomato")) +
  theme_minimal() +
  labs(x = NULL, y = "Count", title = "First Steps participation") +
  theme(legend.position = "none")

p2 <- df %>%
  mutate(low_bwt = factor(low_bwt, labels = c("Normal weight", "Low birth weight"))) %>%
  ggplot(aes(x = low_bwt, fill = low_bwt)) +
  geom_bar(alpha = 0.8, width = 0.4) +
  scale_fill_manual(values = c("steelblue", "tomato")) +
  theme_minimal() +
  labs(x = NULL, y = "Count", title = "Low birth weight outcome") +
  theme(legend.position = "none")

p1 + p2

ggplot(df, aes(x = age)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30,
                 fill = "steelblue", alpha = 0.7) +
  geom_density(color = "navy", linewidth = 0.8) +
  theme_minimal() +
  labs(x = "Maternal age (years)", y = "Density",
       title = "Distribution of maternal age")

df %>%
  mutate(firstep = factor(firstep, labels = c("Non-participant", "Participant")),
         low_bwt = factor(low_bwt, labels = c("Normal", "Low BW"))) %>%
  count(firstep, low_bwt) %>%
  group_by(firstep) %>%
  mutate(prop = n / sum(n)) %>%
  ggplot(aes(x = firstep, y = prop, fill = low_bwt)) +
  geom_col(position = "dodge", alpha = 0.85, width = 0.5) +
  geom_text(aes(label = scales::percent(prop, accuracy = 0.1)),
            position = position_dodge(width = 0.5), vjust = -0.4, size = 3.5) +
  scale_fill_manual(values = c("steelblue", "tomato")) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  theme_minimal() +
  labs(x = NULL, y = "Proportion", fill = NULL,
       title = "Low birth weight rate by First Steps participation")


#############################################
###  Super Learner Setup
#############################################

SL.library <- c("SL.glm", "SL.ranger", "SL.earth", "SL.gam")

# ── Cross-fitting folds ───────────────────────────────────────────────────────
n_folds <- 5
fold_id <- sample(rep(1:n_folds, length.out = N))
folds   <- split(seq_len(N), fold_id)


#############################################
###  Question 1 — AL Linear Regression: First Steps -> Low Birth Weight
#############################################

# ── Nuisance functions ─────────────────────────────────────

y_hat <- numeric(N)
p_hat <- numeric(N)

AL <- data.frame(A, L)

for (k in 1:n_folds) {
  
  test  <- folds[[k]]
  train <- setdiff(seq_len(N), test)
  
  q.model <- SuperLearner(
    Y          = Y[train],
    X          = AL[train, ],
    SL.library = SL.library,
    family     = binomial(),
    method     = "method.NNLS",
    cvControl  = list(V = 5)
  )
  
  p.model <- SuperLearner(
    Y          = A[train],
    X          = data.frame(L[train, ]),  
    SL.library = SL.library,
    family     = binomial(),
    method     = "method.NNLS",
    cvControl  = list(V = 5)
  )
  
  y_hat[test] <- predict(q.model, newdata = AL[test, ])$pred  
  p_hat[test] <- predict(p.model, newdata = data.frame(L[test, ]))$pred
}




# ── Assumption-lean estimator via lm (sandwich SE) ────────────────────────────

# Regress (Y - y_hat) onto (A - p_hat), no intercept
model_q1 <- lm(I(Y - y_hat) ~ -1 + I(A - p_hat))
lm_res_q1 <- coeftest(model_q1, vcov = sandwich)

beta_lm <- lm_res_q1[1, 1]
se_lm   <- lm_res_q1[1, 2]

# ── Assumption-lean estimator via EIF (by hand) ───────────────────────────────

# Point estimate from the efficient influence function
beta_eif <- mean((A - p_hat) * (Y - y_hat)) / mean((A - p_hat)^2)

# Standard error from the EIF
se_eif <- (1 / sqrt(N)) *
  sd((A - p_hat) * (Y - y_hat - beta_eif * (A - p_hat)) / mean((A - p_hat)^2))

# ── Results ───────────────────────────────────────────────────────────────────
results_q1 <- tibble(
  Method    = c("Linear Model", "Eff. Influence Function"),
  Estimate  = c(beta_lm,  beta_eif),
  StdError  = c(se_lm,    se_eif),
  CI_Lower  = c(beta_lm  - 1.96 * se_lm,  beta_eif  - 1.96 * se_eif),
  CI_Upper  = c(beta_lm  + 1.96 * se_lm,  beta_eif  + 1.96 * se_eif)
)

knitr::kable(results_q1, digits = 4,
             caption = "AL linear regression: First Steps on low birth weight")


#############################################
###  Question 2 — Weighted Baseline Characteristics
#############################################

# The naive plug-in estimator suffers plug-in bias: ML estimates p_hat are
# regularized and p*(1-p) is non-linear, so bias does not cancel in the ratio.
# The debiased EIF-based estimator uses (A - p_hat)^2 as weights instead.

L_num <- L   

weighted_summary <- do.call(rbind, lapply(names(L_num), function(varname) {
  X <- L_num[[varname]]
  if (!is.numeric(X)) return(NULL)
  
  w_mean <- mean((A - p_hat)^2 * X) / mean((A - p_hat)^2)
  eif    <- (A - p_hat)^2 * (X - w_mean) / mean(p_hat * (1 - p_hat))
  w_se   <- sd(eif) / sqrt(N)
  
  data.frame(
    variable    = varname,
    w_mean      = round(w_mean, 3),
    w_se        = round(w_se,   3),
    w_CI_lower  = round(w_mean - 1.96 * w_se, 3),
    w_CI_upper  = round(w_mean + 1.96 * w_se, 3)
  )
}))

knitr::kable(weighted_summary,
             col.names = c("Variable", "Weighted mean", "SE", "CI lower", "CI upper"),
             caption   = "Weighted baseline characteristics (overlap weights)")

# ── Classic (unweighted) summary for comparison ───────────────────────────────
classic_summary <- do.call(rbind, lapply(names(L_num), function(varname) {
  X <- L_num[[varname]]
  if (!is.numeric(X)) return(NULL)
  
  c_mean <- mean(X)
  c_se   <- sd(X) / sqrt(N)
  
  data.frame(
    variable    = varname,
    c_mean      = round(c_mean, 3),
    c_se        = round(c_se,   3),
    c_CI_lower  = round(c_mean - 1.96 * c_se, 3),
    c_CI_upper  = round(c_mean + 1.96 * c_se, 3)
  )
}))

knitr::kable(classic_summary,
             col.names = c("Variable", "Classic mean", "SE", "CI lower", "CI upper"),
             caption   = "Classic (unweighted) baseline characteristics")


#############################################
###  Question 3 — AL Linear Regression: Maternal Age -> Low Birth Weight
#############################################

# Exposure is continuous; parity & firstep excluded (potential mediators)

# ── Nuisance functions ─────────────────────────────────────

y2_hat <- numeric(N)
p2_hat <- numeric(N)

al <- data.frame(a, l)

for (k in 1:n_folds) {
  
  test  <- folds[[k]]
  train <- setdiff(seq_len(N), test)
  
  q2.model <- SuperLearner(
    Y          = y[train],
    X          = al[train, ],
    SL.library = SL.library,
    family     = binomial(),
    method     = "method.NNLS",
    cvControl  = list(V = 5)
  )
  
  p2.model <- SuperLearner(
    Y          = a[train],
    X          = data.frame(l[train, ]),  
    SL.library = SL.library,
    family     = gaussian(),
    method     = "method.NNLS",
    cvControl  = list(V = 5)
  )
  
  y2_hat[test] <- predict(q2.model, newdata = al[test, ])$pred  
  p2_hat[test] <- predict(p2.model, newdata = data.frame(l[test, ]))$pred
}

# ── Estimates ─────────────────────────────────────────────────────────────────

# LM with sandwich SE
model_q3  <- lm(I(y - y2_hat) ~ -1 + I(a - p2_hat))
lm_res_q3 <- coeftest(model_q3, vcov = sandwich)

beta_lm3 <- lm_res_q3[1, 1]
se_lm3   <- lm_res_q3[1, 2]

# EIF
beta_eif3 <- mean((a - p2_hat) * (y - y2_hat)) / mean((a - p2_hat)^2)
se_eif3   <- (1 / sqrt(N)) *
  sd((a - p2_hat) * (y - y2_hat - beta_eif3 * (a - p2_hat)) / mean((a - p2_hat)^2))

results_q3 <- tibble(
  Method   = c("Linear Model", "Eff. Influence Function"),
  Estimate = c(beta_lm3, beta_eif3),
  StdError = c(se_lm3,   se_eif3),
  CI_Lower = c(beta_lm3 - 1.96 * se_lm3, beta_eif3 - 1.96 * se_eif3),
  CI_Upper = c(beta_lm3 + 1.96 * se_lm3, beta_eif3 + 1.96 * se_eif3)
)

knitr::kable(results_q3, digits = 6,
             caption = "AL linear regression: maternal age on low birth weight")



#############################################
###  Question 5 — Risk Ratio Scale
#############################################

# Assumption-lean log-linear regression:
# log E(Y^a | L) = alpha(L) + beta*a
# pseudo-outcome: log(Y_hat) - q_log + (Y - Y_hat) / Y_hat

# ── First Steps program ───────────────────────────────────────────────────────

# Step 1: fit E(Y | A, L) using both A and L
AL <- cbind(A = A, L)
Y_hat.mod <- SuperLearner(Y = Y, X = AL,
                          SL.library = SL.library,
                          family = "binomial",
                          cvControl = list(V = 5))
Y_hat <- Y_hat.mod$SL.predict

# Step 2: fit E[log(Y_hat) | L]
log_Y_hat <- log(Y_hat)
q_log.mod <- SuperLearner(Y = log_Y_hat, X = L,
                          SL.library = SL.library,
                          family = "gaussian",
                          cvControl = list(V = 5))
q_log <- q_log.mod$SL.predict

# Step 3: construct pseudo-outcome and regress
pseudo_log <- log_Y_hat - q_log + (Y - Y_hat) / Y_hat
log_mod    <- lm(I(pseudo_log) ~ -1 + I(A - p))
log_res    <- coeftest(log_mod, vcov = sandwich)

log_beta <- log_res[1, 1]
log_se   <- log_res[1, 2]

RR_results <- tibble(
  Method   = "Risk Ratio",
  Estimate = exp(log_beta),
  StdError = log_se,
  CI_Lower = exp(log_beta - 1.96 * log_se),
  CI_Upper = exp(log_beta + 1.96 * log_se)
)

knitr::kable(RR_results, digits = 3,
             caption = "Risk ratio - First Steps program on low birth weight")

# ── Maternal age ──────────────────────────────────────────────────────────────

al <- cbind(a = a, l)
y_hat.mod <- SuperLearner(Y = y, X = al,
                          SL.library = SL.library,
                          family = "binomial",
                          cvControl = list(V = 5))
y_hat <- y_hat.mod$SL.predict

log_y_hat  <- log(y_hat)
q2_log.mod <- SuperLearner(Y = log_y_hat, X = l,
                           SL.library = SL.library,
                           family = "gaussian",
                           cvControl = list(V = 5))
q2_log <- q2_log.mod$SL.predict

pseudo_log2 <- log_y_hat - q2_log + (y - y_hat) / y_hat
log_mod2    <- lm(I(pseudo_log2) ~ -1 + I(a - p2))
log_res2    <- coeftest(log_mod2, vcov = sandwich)

log_beta2 <- log_res2[1, 1]
log_se2   <- log_res2[1, 2]

RR_results2 <- tibble(
  Method   = "Risk Ratio Age",
  Estimate = exp(log_beta2),
  StdError = log_se2,
  CI_Lower = exp(log_beta2 - 1.96 * log_se2),
  CI_Upper = exp(log_beta2 + 1.96 * log_se2)
)

knitr::kable(RR_results2, digits = 4,
             caption = "Risk ratio - maternal age on low birth weight")


#############################################
###  Question 6 — Odds Ratio Scale
#############################################

# Assumption-lean logistic regression:
# logit E(Y^a | L) = alpha(L) + beta*a
# pseudo-outcome: logit(Y_hat) - q_logit + (Y - Y_hat) / (Y_hat*(1-Y_hat))

# ── First Steps program ───────────────────────────────────────────────────────

logit_Y_hat <- log(Y_hat / (1 - Y_hat))
q_logit.mod <- SuperLearner(Y = logit_Y_hat, X = L,
                            SL.library = SL.library,
                            family = "gaussian",
                            cvControl = list(V = 5))
q_logit <- q_logit.mod$SL.predict

pseudo_logit <- logit_Y_hat - q_logit + (Y - Y_hat) / (Y_hat * (1 - Y_hat))
logit_mod    <- lm(I(pseudo_logit) ~ -1 + I(A - p))
logit_res    <- coeftest(logit_mod, vcov = sandwich)

logit_beta <- logit_res[1, 1]
logit_se   <- logit_res[1, 2]

OR_results <- tibble(
  Method   = "Odds Ratio",
  Estimate = exp(logit_beta),
  StdError = logit_se,
  CI_Lower = exp(logit_beta - 1.96 * logit_se),
  CI_Upper = exp(logit_beta + 1.96 * logit_se)
)

knitr::kable(OR_results, digits = 3,
             caption = "Odds ratio - First Steps program on low birth weight")

# ── Maternal age ──────────────────────────────────────────────────────────────

logit_y_hat  <- log(y_hat / (1 - y_hat))
q2_logit.mod <- SuperLearner(Y = logit_y_hat, X = l,
                             SL.library = SL.library,
                             family = "gaussian",
                             cvControl = list(V = 5))
q2_logit <- q2_logit.mod$SL.predict

pseudo_logit2 <- logit_y_hat - q2_logit + (y - y_hat) / (y_hat * (1 - y_hat))
logit_mod2    <- lm(I(pseudo_logit2) ~ -1 + I(a - p2))
logit_res2    <- coeftest(logit_mod2, vcov = sandwich)

logit_beta2 <- logit_res2[1, 1]
logit_se2   <- logit_res2[1, 2]

OR_results2 <- tibble(
  Method   = "Odds Ratio Age",
  Estimate = exp(logit_beta2),
  StdError = logit_se2,
  CI_Lower = exp(logit_beta2 - 1.96 * logit_se2),
  CI_Upper = exp(logit_beta2 + 1.96 * logit_se2)
)

knitr::kable(OR_results2, digits = 3,
             caption = "Odds ratio - maternal age on low birth weight")

# ── Visualizing the odds-ratio transformation ─────────────────────────────────
expit <- function(x) 1 / (1 + exp(-x))
logit <- function(p) log(p / (1 - p))

x_seq <- seq(0.001, 0.999, length.out = 1000)
y_seq <- expit(logit(x_seq) + logit_beta)

plot(x_seq, y_seq, type = "l", col = "steelblue", lwd = 2,
     xlab = expression(E(Y ~ "|" ~ A==0 ~ "," ~ L)),
     ylab = expression(E(Y ~ "|" ~ A==1 ~ "," ~ L)),
     main = expression("Implied" ~ E(Y ~ "|" ~ A==1 ~ "," ~ L) ~
                         "= expit(logit(" * hat(Y)[A==0] * ") + " * hat(beta) * ")"))
abline(0, 1, col = "gray50", lty = 2)
legend("topleft", legend = c("OR model", "Identity"),
       col = c("steelblue", "gray50"), lty = c(1, 2), bty = "n")


#############################################
###  Summary of Results
#############################################

summary_all <- tibble(
  Exposure = rep(c("First Steps", "Maternal Age"), each = 3),
  Scale    = rep(c("Risk Difference", "Risk Ratio", "Odds Ratio"), times = 2),
  Estimate = c(beta_eif,        exp(log_beta),   exp(logit_beta),
               beta_eif3,       exp(log_beta2),  exp(logit_beta2)),
  SE       = c(se_eif,          log_se,          logit_se,
               se_eif3,         log_se2,         logit_se2),
  Null     = rep(c(0, 1, 1), 2)
)

summary_all <- summary_all %>%
  mutate(CI_Lower = if_else(Scale == "Risk Difference",
                            Estimate - 1.96 * SE, exp(log(Estimate) - 1.96 * SE)),
         CI_Upper = if_else(Scale == "Risk Difference",
                            Estimate + 1.96 * SE, exp(log(Estimate) + 1.96 * SE)))

knitr::kable(summary_all %>%
               select(-Null) %>%
               mutate(across(where(is.numeric), ~round(.x, 4))),
             caption = "Summary of all assumption-lean regression results")

# ── Forest plot ───────────────────────────────────────────────────────────────
summary_all %>%
  mutate(Label = paste(Exposure, "-", Scale),
         Label = factor(Label, levels = rev(unique(Label)))) %>%
  ggplot(aes(x = Estimate, y = Label, color = Exposure)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0.3, linewidth = 0.8) +
  geom_vline(aes(xintercept = Null), linetype = "dashed", color = "gray50") +
  facet_wrap(~ Scale, scales = "free_x", nrow = 1) +
  scale_color_manual(values = c("steelblue", "tomato")) +
  theme_minimal() +
  labs(x = "Estimate (95% CI)", y = NULL, color = "Exposure",
       title = "Assumption-Lean Regression - King County 2001",
       subtitle = "Risk difference, risk ratio, and odds ratio scales") +
  theme(legend.position = "bottom")