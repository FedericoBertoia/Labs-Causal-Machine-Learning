

Supplementary: Code

### 0 Setup

library(hdm)           # pension data
library(SuperLearner)  # ensemble nuisance estimation
library(tmle)          # TMLE implementation
library(ggplot2)       # plotting
library(sandwich)      # robust SEs

remotes::install_github("ehkennedy/npcausal")
library(npcausal)

#Load some functions from the Lab 5 solutions:

# Propensity-score cross-fitting 
estimate_pi_nonpar <- function(A, X, sl_lib, nsplits = 5, family = binomial()) {
  N <- length(A)
  folds <- sample(rep(seq_len(nsplits), length.out = N))
  preds <- numeric(N)
  coef_sum <- rep(0, length(sl_lib))
  for (v in seq_len(nsplits)) {
    train <- which(folds != v)
    test  <- which(folds == v)
    fit <- SuperLearner(
      Y           = A[train],
      X           = X[train, , drop=FALSE],
      newX        = X[test,  , drop=FALSE],
      SL.library  = sl_lib,
      family      = family,
      method      = "method.NNloglik"
    )
    preds[test]    <- fit$SL.predict
    coef_sum       <- coef_sum + fit$coef
  }
  list(preds = preds, coefs = coef_sum / nsplits)
}

estimate_mu_nonpar <- function(Y, A, X, sl_lib, nsplits = 5, family = gaussian()) {
  N <- length(Y)
  folds <- sample(rep(seq_len(nsplits), length.out = N))
  mu0 <- numeric(N); mu1 <- numeric(N)
  coef0 <- coef1 <- rep(0, length(sl_lib))
  for (v in seq_len(nsplits)) {
    train <- which(folds != v)
    test  <- which(folds == v)
    # Fit on treated = 1
    fit1 <- SuperLearner(
      Y           = Y[A[train]==1 & train],
      X           = X[A[train]==1 & train, , drop=FALSE],
      newX        = X[test, , drop=FALSE],
      SL.library  = sl_lib,
      family      = family
    )
    mu1[test] <- fit1$SL.predict
    coef1     <- coef1 + fit1$coef
    # Fit on treated = 0
    fit0 <- SuperLearner(
      Y           = Y[A[train]==0 & train],
      X           = X[A[train]==0 & train, , drop=FALSE],
      newX        = X[test, , drop=FALSE],
      SL.library  = sl_lib,
      family      = family
    )
    mu0[test] <- fit0$SL.predict
    coef0     <- coef0 + fit0$coef
  }
  list(mu0  = mu0,
       mu1  = mu1,
       coefs = list(coef0 / nsplits, coef1 / nsplits))
}

### 1 Load data and define objects

data(pension)

Y      <- pension$net_tfa
A      <- pension$e401
L      <- pension[, c("age","inc","fsize","educ","marr","twoearn","db","ira","hown")]
sl_lib <- c("SL.glmnet", "SL.gam")

### 2 Unadjusted ATE 

unadj_ATE <- mean(Y[A==1]) - mean(Y[A==0])
unadj_SE  <- sqrt(var(Y[A==1])/sum(A==1) + var(Y[A==0])/sum(A==0))
cat("Unadjusted ATE =", round(unadj_ATE,1), "SE =", round(unadj_SE,1), "\n")

### 3. AIPW via npcausal::ate()

aipw_res <- ate(
  y       = Y,
  a       = A,
  x       = L,
  nsplits = 5,
  sl.lib  = sl_lib,
  trim    = NULL
)
print(aipw_res$res)  # E[Y(0)], E[Y(1)], ATE  SE & 95% CI
cat("\n")

### 4 Manual cross-fitting + AIPW 

# 4a) Propensity scores g(L)
ghat <- estimate_pi_nonpar(A, L, sl_lib, nsplits=5)$preds

# 4b) Outcome regressions Q(L), Q(L)
mu_fit <- estimate_mu_nonpar(Y, A, L, sl_lib, nsplits=5)
Q0     <- mu_fit$mu0
Q1     <- mu_fit$mu1
tau    <- Q1 - Q0

# 4c) Influence-curve and ATE
N      <- length(Y)
IC     <- (A/ghat - (1-A)/(1-ghat)) * (Y - ifelse(A==1, Q1, Q0)) + tau
ate_m  <- mean(IC)
se_m   <- sd(IC)/sqrt(N)
cat("Manual AIPW  ATE =", round(ate_m,1), "SE =", round(se_m,1), "\n\n")

### 5 Diagnostics for AIPW

# 5a. Density of g(L) by A
df_propensity_scores <- data.frame(
  ghat = ghat,
  A    = factor(A, levels = c(0,1), labels = c("Ineligible","Eligible"))
)
ggplot(df_propensity_scores, aes(x = ghat, color = A)) +
  geom_density() +
  labs(
    title = "Density of Propensity Scores g(L) by Eligibility",
    x     = expression(hat(g)(L)),
    color = "Eligibility"
  )

# 5b. Boxplot of inverse-propensity weights
df_weights <- data.frame(
  w = ifelse(A == 1, 1/ghat, 1/(1 - ghat)),
  A = factor(A, levels = c(0,1), labels = c("Ineligible","Eligible"))
)
ggplot(df_weights, aes(x = A, y = w)) +
  geom_boxplot() +
  labs(
    title = "Inverse-Propensity Weights by Eligibility",
    x     = "Eligibility",
    y     = expression(1 / hat(g)(L) ~ "or" ~ 1/(1 - hat(g)(L)))
  )

summary(IC)

### 6 Covariates and weights

# 1. Build DF with covariates, propensity & weight
df <- data.frame(
  L,
  ghat   = ghat,
  weight = ifelse(A == 1, 1/ghat, 1/(1 - ghat))
)

# 2. Restrict to the eligible (A=1) households
df_e <- subset(df, A == 1)

# 3a. See the top 5 highestweight treated units
top5 <- df_e[order(-df_e$weight), ][1:5, ]
print(top5)

# 3b. Quick correlations between each L_j and the weight
cors <- sapply(names(L), function(var){
  cor(df_e[[var]], df_e$weight)
})
sort(cors, decreasing=TRUE)

### 7 AIPW ATT via npcausal::att() 

att_res <- att(
  y       = Y,
  a       = A,
  x       = L,
  nsplits = 5,
  sl.lib  = sl_lib)

### 8 TMLE ATE

tmle_fit <- tmle(
  Y               = Y,
  A               = A,
  W               = L,
  Q.SL.library    = sl_lib,
  g.SL.library    = sl_lib,
  family          = "gaussian",
  V.Q             = 5,
  V.g             = 5
)

# Extract and print ATE
ate_tmle    <- tmle_fit$estimates$ATE$psi
se_ate_tmle <- sqrt(tmle_fit$estimates$ATE$var.psi)

### 9 TMLE ATT

# 1) Obtain fitted g(L) and Q0(L)
ghat  <- estimate_pi_nonpar(A, L, sl_lib, nsplits = 5)$preds
mu0_1 <- estimate_mu_nonpar(Y, A, L, sl_lib, nsplits = 5)$mu0

# 2) Fluctuation model on the control group
fluct_mod <- lm(
  Y[A == 0] ~ -1 +
    offset(mu0_1[A == 0]) +
    I(ghat[A == 0] / (1 - ghat[A == 0]))
)
epsilon <- coef(fluct_mod)

# 3) Update Q1 via targeted fluctuation
Q1_star <- mu0_1 + epsilon * (ghat / (1 - ghat))

# 4) Compute ATT and its SE via the influence curve
pA       <- mean(A)
ATT_tmle <- mean((A / pA) * (Y - Q1_star))
IC_att   <- (A / pA) * (Y - Q1_star - ATT_tmle) -
  ((1 - A) * ghat / (1 - ghat)) / pA * (Y - Q1_star)
se_att_tmle <- sd(IC_att) / sqrt(length(Y))

cat(sprintf("TMLE ATT = %.1f (SE = %.1f)\n",
            ATT_tmle, se_att_tmle))

### 9 E{Cov(A,Y|L)}=0

# 7a) Crossfit E[Y | L] with SuperLearner

#Copy-paste lab 5 but without A 
estimate_m_nonpar <- function(Y, X, sl_lib, nsplits = 5, family = gaussian()) {
  N <- length(Y)
  folds <- sample(rep(seq_len(nsplits), length.out = N))
  preds <- numeric(N)
  coef_sum <- rep(0, length(sl_lib))
  for (v in seq_len(nsplits)) {
    train <- which(folds != v)
    test  <- which(folds == v)
    fit <- SuperLearner(
      Y          = Y[train],
      X          = X[train, , drop = FALSE],
      newX       = X[test,  , drop = FALSE],
      SL.library = sl_lib,
      family     = family
    )
    preds[test]  <- fit$SL.predict
    coef_sum     <- coef_sum + fit$coef
  }
  list(preds = preds, coefs = coef_sum / nsplits)
}

# (re)compute g(L) if needed
g_out <- estimate_pi_nonpar(A, L, sl_lib, nsplits = 5)
ghat  <- g_out$preds

# compute m(L) = E[Y | L]
m_out <- estimate_m_nonpar(Y, L, sl_lib, nsplits = 5)
mhat  <- m_out$preds

# form Di = (A - g(L)) * (Y - m(L))
D     <- (A - ghat) * (Y - mhat)
N     <- length(D)

# estimate  and its SE
theta_hat <- mean(D)
se_hat    <- sd(D) / sqrt(N)

# twosided pvalue for H0:  = 0
z_stat <- theta_hat / se_hat
pval   <- 2 * (1 - pnorm(abs(z_stat)))





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Google Colab R Project Script for Causal Machine Learning
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This script includes package installation for a Colab environment.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# --- 0.A Package Installation for Colab Environment ---
# This section will install necessary packages. It might take a few minutes.

if (!requireNamespace("hdm", quietly = TRUE)) install.packages("hdm")
if (!requireNamespace("SuperLearner", quietly = TRUE)) install.packages("SuperLearner")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("gbm", quietly = TRUE)) install.packages("gbm")
if (!requireNamespace("tmle", quietly = TRUE)) install.packages("tmle")
if (!requireNamespace("Metrics", quietly = TRUE)) install.packages("Metrics") # For optional MSE check in Q5
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("randomForest", quietly = TRUE)) install.packages("randomForest")

# rlearner might need to be installed from GitHub
if (!requireNamespace("rlearner", quietly = TRUE)) {
  devtools::install_github("xnie/rlearner")
}

# rlearner might need to be installed from GitHub
if (!requireNamespace("rlearner", quietly = TRUE)) {
  devtools::install_github("xnie/rlearner")
}

# --- 0.B Setup ---
set.seed(12345) # For reproducibility

library(hdm)
library(SuperLearner)
library(dplyr)
library(ggplot2)
library(rlearner)
library(gbm)
library(tmle)
library(Metrics)

# Load data
data(pension) # Full dataset: 9915 observations

# Define global parameters
K_FOLDS <- 5 # Number of folds for cross-fitting

outcome_var <- "net_tfa"
treatment_var <- "e401" # Eligibility
covariate_names <- c("age", "inc", "fsize", "educ", "marr", "twoearn", "db", "pira", "hown")

X_df <- pension[, covariate_names]
Y_vec <- pension[[outcome_var]]
A_vec <- pension[[treatment_var]]

SL_LIBRARY_GAUSSIAN <- c("SL.mean", "SL.glm", "SL.glmnet", "SL.randomForest", "SL.gbm")
SL_LIBRARY_BINOMIAL <- c("SL.mean", "SL.glm", "SL.glmnet", "SL.randomForest", "SL.gbm")

bound_pi <- function(pi_hat, bounds = c(0.025, 0.975)) {
  pi_hat_bounded <- pmin(pmax(pi_hat, bounds[1]), bounds[2])
  return(pi_hat_bounded)
}

# --- Question 2: Estimate Eligibility Effect (AIPW & TMLE) ---
cat("--- Question 2: AIPW & TMLE --- \n")

unadjusted_est_q2 <- mean(Y_vec[A_vec == 1]) - mean(Y_vec[A_vec == 0])
unadjusted_ttest_q2 <- t.test(Y_vec ~ A_vec)
cat("Unadjusted ATE (e401 on net_tfa):\n")
print(unadjusted_ttest_q2)

folds <- sample(rep(1:K_FOLDS, length.out = nrow(pension)))
mu_hat_0_cv <- numeric(nrow(pension))
mu_hat_1_cv <- numeric(nrow(pension))
mu_hat_A_cv <- numeric(nrow(pension))
pi_hat_cv <- numeric(nrow(pension))

cat("Starting nuisance function estimation for Q2 (this may take time)...\n")
for (k in 1:K_FOLDS) {
  train_idx <- which(folds != k)
  test_idx <- which(folds == k)
  X_train <- X_df[train_idx, ]; Y_train <- Y_vec[train_idx]; A_train <- A_vec[train_idx]
  X_test <- X_df[test_idx, ]

  pi_sl_fit <- SuperLearner(Y = A_train, X = X_train, family = binomial(), SL.library = SL_LIBRARY_BINOMIAL, cvControl = list(V = K_FOLDS))
  pi_hat_cv[test_idx] <- predict(pi_sl_fit, newdata = X_test)$pred[,1]

  mu1_sl_fit <- SuperLearner(Y = Y_train[A_train == 1], X = X_train[A_train == 1, ], family = gaussian(), SL.library = SL_LIBRARY_GAUSSIAN, cvControl = list(V = K_FOLDS))
  mu_hat_1_cv[test_idx] <- predict(mu1_sl_fit, newdata = X_test)$pred[,1]

  mu0_sl_fit <- SuperLearner(Y = Y_train[A_train == 0], X = X_train[A_train == 0, ], family = gaussian(), SL.library = SL_LIBRARY_GAUSSIAN, cvControl = list(V = K_FOLDS))
  mu_hat_0_cv[test_idx] <- predict(mu0_sl_fit, newdata = X_test)$pred[,1]

  mu_hat_A_cv[test_idx] <- ifelse(A_vec[test_idx] == 1, mu_hat_1_cv[test_idx], mu_hat_0_cv[test_idx])
  cat("Completed fold", k, "of", K_FOLDS, "for Q2 nuisance functions.\n")
}
cat("Nuisance function estimation for Q2 complete.\n")

pi_hat_cv_bounded <- bound_pi(pi_hat_cv)

term1_aipw <- (A_vec / pi_hat_cv_bounded) * (Y_vec - mu_hat_1_cv)
term2_aipw <- ((1 - A_vec) / (1 - pi_hat_cv_bounded)) * (Y_vec - mu_hat_0_cv)
aipw_influence_curve <- term1_aipw - term2_aipw + mu_hat_1_cv - mu_hat_0_cv
ate_aipw <- mean(aipw_influence_curve)
se_aipw <- sd(aipw_influence_curve) / sqrt(nrow(pension))
ci_aipw <- c(ate_aipw - 1.96 * se_aipw, ate_aipw + 1.96 * se_aipw)
cat("AIPW ATE (e401 on net_tfa):\n")
cat("Estimate:", ate_aipw, "SE:", se_aipw, "95% CI: [", ci_aipw[1], ",", ci_aipw[2], "]\n\n")

Q_init <- cbind(mu_hat_0_cv, mu_hat_1_cv)
g_init <- pi_hat_cv_bounded
W_tmle <- X_df
tmle_fit_q2 <- tmle(Y = Y_vec, A = A_vec, W = W_tmle, Q = Q_init, g1W = g_init, family = "gaussian")
cat("TMLE ATE (e401 on net_tfa):\n")
print(tmle_fit_q2)
cat("\n")

# --- Question 3: Test E{Cov(A,Y|L)}=0 ---
cat("--- Question 3: Test E{Cov(A,Y|L)}=0 --- \n")
mu_X_hat_cv <- numeric(nrow(pension))
cat("Starting E[Y|X] estimation for Q3 (this may take time)...\n")
for (k in 1:K_FOLDS) {
  train_idx <- which(folds != k); test_idx <- which(folds == k)
  X_train <- X_df[train_idx, ]; Y_train <- Y_vec[train_idx]; X_test <- X_df[test_idx, ]
  mu_X_sl_fit <- SuperLearner(Y = Y_train, X = X_train, family = gaussian(), SL.library = SL_LIBRARY_GAUSSIAN, cvControl = list(V = K_FOLDS))
  mu_X_hat_cv[test_idx] <- predict(mu_X_sl_fit, newdata = X_test)$pred[,1]
  cat("Completed fold", k, "of", K_FOLDS, "for Q3 E[Y|X].\n")
}
cat("E[Y|X] estimation for Q3 complete.\n")

resid_A_q3 <- A_vec - pi_hat_cv_bounded
resid_Y_q3 <- Y_vec - mu_X_hat_cv
theta_q3_est <- mean(resid_A_q3 * resid_Y_q3)
phi_q3 <- resid_A_q3 * resid_Y_q3 - theta_q3_est
var_phi_q3 <- var(phi_q3)
se_theta_q3 <- sqrt(var_phi_q3 / nrow(pension))
z_score_q3 <- theta_q3_est / se_theta_q3
p_value_q3 <- 2 * pnorm(-abs(z_score_q3))
cat("Estimate of E{Cov(A,Y|L)} (theta):", theta_q3_est, "\n")
cat("SE(theta):", se_theta_q3, "\n")
cat("Z-score:", z_score_q3, "P-value:", p_value_q3, "\n\n")

# --- Question 4: R-learner for Heterogeneous Effects ---
cat("--- Question 4: R-learner --- \n")
cat("Starting R-learner estimation (this may take time)...\n")
X_mat <- as.matrix(X_df)
rboost_fit_q4 <- rboost(X_mat, A_vec, Y_vec)
cate_estimates_q4 <- predict(rboost_fit_q4, X_mat)
cat("R-learner estimation complete.\n")

pension_q4_results <- pension
pension_q4_results$cate_rlearner <- cate_estimates_q4

hist_cate_q4 <- ggplot(pension_q4_results, aes(x = cate_rlearner)) +
  geom_histogram(bins = 50, fill = "skyblue", color = "black", aes(y=..density..)) +
  geom_density(alpha=0.5, fill="lightblue") +
  labs(title = "Histogram of Estimated CATEs (R-learner)", x = "Estimated CATE ((X))", y = "Density") + theme_minimal()
print(hist_cate_q4)

plot_cate_age <- ggplot(pension_q4_results, aes(x = age, y = cate_rlearner)) +
  geom_point(alpha = 0.3) + geom_smooth(method = "loess", se = FALSE, color = "blue") +
  labs(title = "CATE vs. Age", x = "Age", y = "Estimated CATE ((X))") + theme_minimal()
print(plot_cate_age)

plot_cate_inc <- ggplot(pension_q4_results, aes(x = inc, y = cate_rlearner)) +
  geom_point(alpha = 0.3) + geom_smooth(method = "loess", se = FALSE, color = "blue") +
  labs(title = "CATE vs. Income", x = "Income", y = "Estimated CATE ((X))") + theme_minimal()
print(plot_cate_inc)

plot_cate_educ <- ggplot(pension_q4_results, aes(x = educ, y = cate_rlearner)) +
  geom_point(alpha = 0.3) + geom_smooth(method = "loess", se = FALSE, color = "blue") +
  labs(title = "CATE vs. Education", x = "Years of Education", y = "Estimated CATE ((X))") + theme_minimal()
print(plot_cate_educ)
cat("R-learner CATE plots generated.\n\n")

# --- Question 5: Predicting Y and Estimating Prediction Error ---
cat("--- Question 5: Predicting Y^1 and Error --- \n")
Ci_q5 <- mu_hat_1_cv + (A_vec / pi_hat_cv_bounded) * (Y_vec - mu_hat_1_cv)
X_small_df <- pension[, c("age", "inc", "educ")]

Ytilde1_cv <- numeric(nrow(pension))
cat("Starting Ytilde1 DR-learner prediction model estimation (this may take time)...\n")
for (k in 1:K_FOLDS) {
  train_idx <- which(folds != k); test_idx <- which(folds == k)
  Ci_train <- Ci_q5[train_idx]; X_small_train <- X_small_df[train_idx, ]; X_small_test <- X_small_df[test_idx, ]

  sl_ytilde1_fit <- SuperLearner(Y = Ci_train, X = X_small_train, family = gaussian(), SL.library = SL_LIBRARY_GAUSSIAN, cvControl = list(V = K_FOLDS))
  Ytilde1_cv[test_idx] <- predict(sl_ytilde1_fit, newdata = X_small_test)$pred[,1]
  cat("Completed fold", k, "of", K_FOLDS, "for Q5 Ytilde1 model.\n")
}
cat("Ytilde1 DR-learner prediction model estimation complete.\n")

pension_q5_results <- pension
pension_q5_results$Ytilde1 <- Ytilde1_cv
hist_Ytilde1_q5 <- ggplot(pension_q5_results, aes(x = Ytilde1)) +
  geom_histogram(bins = 50, fill = "coral", color = "black", aes(y=..density..)) +
  geom_density(alpha=0.5, fill="orange") +
  labs(title = "Histogram of Predicted Y^1 (Y)", x = "Predicted Y^1", y = "Density") + theme_minimal()
print(hist_Ytilde1_q5)

cat("Starting E[(Y1-Ytilde1)^2] error estimation (this may take time)...\n")
D_values_q5_all <- (Y_vec - Ytilde1_cv)^2 # (Y_i - Y_i)^2 for all i
mu_D_hat_cv <- numeric(nrow(pension))

for (k in 1:K_FOLDS) {
  train_idx_fold_k <- which(folds != k)
  A1_train_fold_k_idx_original <- intersect(which(A_vec == 1), train_idx_fold_k)

  D_train_fold_k <- D_values_q5_all[A1_train_fold_k_idx_original] # Use (Y-Ytilde)^2 from A=1 in train
  X_full_A1_train_fold_k <- X_df[A1_train_fold_k_idx_original, ]

  test_idx_fold_k <- which(folds == k)
  X_full_test_fold_k <- X_df[test_idx_fold_k, ]

  if(nrow(X_full_A1_train_fold_k) > 1 && length(D_train_fold_k) > 1 && var(D_train_fold_k, na.rm=TRUE) > 1e-6 ) { # Added var check
      sl_mu_D_fit <- SuperLearner(Y = D_train_fold_k, X = X_full_A1_train_fold_k, newX = X_full_test_fold_k, family = gaussian(), SL.library = SL_LIBRARY_GAUSSIAN, cvControl = list(V = K_FOLDS))
      mu_D_hat_cv[test_idx_fold_k] <- predict(sl_mu_D_fit, newdata = X_full_test_fold_k)$pred[,1]
  } else {
      # Fallback if too few A=1 in a fold, no variance in D, or SL fails
      # Global mean of D only among treated
      D_global_mean_treated <- mean( (Y_vec[A_vec==1] - Ytilde1_cv[A_vec==1])^2 , na.rm=TRUE)
      mu_D_hat_cv[test_idx_fold_k] <- D_global_mean_treated
  }
   cat("Completed fold", k, "of", K_FOLDS, "for Q5 mu_D model.\n")
}
cat("E[(Y1-Ytilde1)^2] error component mu_D estimation complete.\n")

term_A_div_pi <- ifelse(pi_hat_cv_bounded > 1e-6 & A_vec == 1, A_vec / pi_hat_cv_bounded, 0)
if_term_error_q5 <- term_A_div_pi * (D_values_q5_all - mu_D_hat_cv) + mu_D_hat_cv
theta_err_est_q5 <- mean(if_term_error_q5)
se_theta_err_q5 <- sd(if_term_error_q5 - theta_err_est_q5) / sqrt(nrow(pension))
ci_theta_err_q5 <- c(theta_err_est_q5 - 1.96 * se_theta_err_q5, theta_err_est_q5 + 1.96 * se_theta_err_q5)

cat("Estimated Counterfactual Prediction Error E{(Y^1 - Y)^2}:\n")
cat("Estimate:", theta_err_est_q5, "SE:", se_theta_err_q5, "95% CI: [", ci_theta_err_q5[1], ",", ci_theta_err_q5[2], "]\n\n")

# --- Final Remarks ---
cat("--- End of Analysis --- \n")
hist(pi_hat_cv_bounded, main='Distribution of Propensity Scores P(A=1|L)', xlab="P(A=1|L)")