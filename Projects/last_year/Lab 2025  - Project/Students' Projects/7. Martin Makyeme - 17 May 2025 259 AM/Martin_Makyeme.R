# Load necessary libraries
library(hdm)             # for the pension dataset
library(SuperLearner)    # For flexible ML-based regression
library(tmle)            # for TMLE and AIPW estimation
library(ggplot2)        # For plotting
library(patchwork)      # for combining multiple plots
library(boot)        # # For estimating confidence intervals using bootstrap resampling

############ Task 2 #############
set.seed(123)  # to ensure reproducibility

#Loading the pension dataset
data(pension)

# checking coloumn names
names(pension)

# Keeping only relevant variables for the analysis
vars <- c("net_tfa", "e401", "p401", "age", "inc", "fsize", "educ", 
          "male", "marr", "twoearn", "db", "pira", "hown")
data <- na.omit(pension[, vars])  # Drop rows with missing values

#Define outcome (Y), treatment (A), and covariates (X)
Y <- data$net_tfa           # Outcome variable: net financial assets
A <- data$e401              # treatment variable: eligibility for 401(k)
X <- data[, !(names(data) %in% c("net_tfa", "e401", "p401"))]  # Covariates

#Specifying learners to be used in SuperLearner
SL.library <- c("SL.glm", "SL.glmnet", "SL.ranger")  # Include linear, penalized, and tree-based learners

# Run TMLE (which also includes AIPW estimates internally)
tmle.out <- tmle(Y = Y,
                 A = A,
                 W = X,
                 family = "gaussian",            # continuous outcome
                 Q.SL.library = SL.library,      # super learner for outcome regression
                 g.SL.library = SL.library)      # Super learner for propensity score

# printing TMLE estimate of ATE and its 95% confidence interval
cat("ATE Estimate (TMLE):", tmle.out$estimates$ATE$psi, "\\n")
cat("95% CI (TMLE): [", tmle.out$estimates$ATE$CI[1], ", ", tmle.out$estimates$ATE$CI[2], "]\\n")

#Printing AIPW estimate and its confidence interval
cat("ATE Estimate (AIPW):", tmle.out$estimates$ATE$DoublyRobust$est, "\\n")
cat("95% CI (AIPW): [", tmle.out$estimates$ATE$DoublyRobust$CI[1], ", ", tmle.out$estimates$ATE$DoublyRobust$CI[2], "]\\n")

# showing naive difference in means without adjustment
unadjusted <- mean(Y[A == 1]) - mean(Y[A == 0])
cat("Unadjusted Difference in Means:", unadjusted, "\\n")


############# Task 3 #######################################################

# Define variables
vars <- c("net_tfa", "e401", "age", "inc", "fsize", "educ",
          "male", "marr", "twoearn", "db", "pira", "hown")
data <- na.omit(pension[, vars])
Y <- data$net_tfa
A <- data$e401
L <- data[, !(names(data) %in% c("net_tfa", "e401"))]

n <- nrow(data)
K <- 5
folds <- sample(rep(1:K, length.out = n))

e_hat <- numeric(n)
mu_hat <- numeric(n)

SL.library <- c("SL.glm", "SL.glmnet", "SL.ranger")

#Cross-fitting loop
for (k in 1:K) {
  train_idx <- which(folds != k)
  test_idx <- which(folds == k)
  
  # estimate E[A | L] on training fold
  ps_model <- SuperLearner(Y = A[train_idx], X = L[train_idx, ], family = binomial(), SL.library = SL.library)
  e_hat[test_idx] <- predict(ps_model, newdata = L[test_idx, ])$pred
  
  # Estimate E[Y | L] on training fold
  mu_model <- SuperLearner(Y = Y[train_idx], X = L[train_idx, ], family = gaussian(), SL.library = SL.library)
  mu_hat[test_idx] <- predict(mu_model, newdata = L[test_idx, ])$pred
}

#Estimate theta
theta_hat <- mean((A - e_hat) * (Y - mu_hat))

# Estimate EIF
phi <- (A - e_hat) * (Y - mu_hat) - theta_hat
se_theta <- sd(phi) / sqrt(n)
z_score <- theta_hat / se_theta
p_value <- 2 * (1 - pnorm(abs(z_score)))

# Output
cat("Estimated θ = E[Cov(A, Y | L)]:", round(theta_hat, 2), "\n")
cat("Standard Error:", round(se_theta, 2), "\n")
cat("Z-score:", round(z_score, 2), "\n")
cat("P-value for H0: θ = 0:", signif(p_value, 3), "\n")

hist(phi, breaks = 50, main = "Distribution of EIF Values", xlab = "Influence Function φ", col = "lightblue")

####### Task 4 ##############################################################################################

# Load and clean the data
data(pension)
df <- pension
df <- df[, c("net_tfa", "e401", "age", "inc", "educ", 
             "fsize", "marr", "twoearn", "db", "pira", "hown")]
df <- df[complete.cases(df), ]  # Drop rows with missing values

#Define variables
Y <- df$net_tfa              # Outcome: net total financial assets
A <- df$e401                 # Treatment: 401(k) eligibility
X <- df[, c("age", "inc", "educ", "fsize", "marr", 
            "twoearn", "db", "pira", "hown")]  # Covariates

# define stable SuperLearner library to reduce errors
SL.lib <- c("SL.mean", "SL.glm")

# Wrapper to safely handle SuperLearner predictions
safe_predict_SL <- function(model, newdata) {
  tryCatch({
    p <- predict(model, newdata = newdata, onlySL = TRUE)
    if (is.null(p$pred)) stop("Prediction failed")
    as.vector(p$pred)
  }, error = function(e) {
    warning("SuperLearner prediction failed. Returning NA.")
    rep(NA, nrow(newdata))
  })
}

#set up 5-fold cross-fitting
set.seed(123)
K <- 5
folds <- sample(rep(1:K, length.out = nrow(df)))
cate_estimates <- rep(NA, nrow(df))  # To store estimated CATEs

# Cross-fitting loop
for (fold in 1:K) {
  idx_test <- which(folds == fold)
  idx_train <- which(folds != fold)
  
  # Split training data into two halves
  idx_A <- idx_train[1:floor(length(idx_train)/2)]
  idx_B <- idx_train[(floor(length(idx_train)/2)+1):length(idx_train)]
  
  # Train on A, predict on B
  mu_A <- SuperLearner(Y = Y[idx_A], X = X[idx_A,], SL.library = SL.lib)
  ps_A <- SuperLearner(Y = A[idx_A], X = X[idx_A,], family = binomial(), SL.library = SL.lib)
  mu_B <- safe_predict_SL(mu_A, X[idx_B,])
  ps_B <- safe_predict_SL(ps_A, X[idx_B,])
  Y_resid_B <- Y[idx_B] - mu_B
  A_resid_B <- A[idx_B] - ps_B
  
  # Train on B, predict on A
  mu_Bmod <- SuperLearner(Y = Y[idx_B], X = X[idx_B,], SL.library = SL.lib)
  ps_Bmod <- SuperLearner(Y = A[idx_B], X = X[idx_B,], family = binomial(), SL.library = SL.lib)
  mu_A_pred <- safe_predict_SL(mu_Bmod, X[idx_A,])
  ps_A_pred <- safe_predict_SL(ps_Bmod, X[idx_A,])
  Y_resid_A <- Y[idx_A] - mu_A_pred
  A_resid_A <- A[idx_A] - ps_A_pred
  
  #construct pseudo-outcomes and weights
  pseudo_Y <- c(Y_resid_A / (A_resid_A + 1e-8), Y_resid_B / (A_resid_B + 1e-8))
  weights <- c(A_resid_A^2, A_resid_B^2)
  X_combined <- rbind(X[idx_A,], X[idx_B,])
  
  #Estimate CATE using weighted regression if valid
  valid <- is.finite(pseudo_Y) & weights > 1e-6
  if (sum(valid) > 30) {
    r_model <- SuperLearner(Y = pseudo_Y[valid], X = X_combined[valid,],
                            SL.library = SL.lib, obsWeights = weights[valid])
    cate_estimates[idx_test] <- safe_predict_SL(r_model, X[idx_test,])
  } else {
    warning(paste("Fold", fold, "skipped: insufficient valid samples"))
  }
}

# Plot distribution of CATEs
hist(cate_estimates,
     main = "Distribution of Individual Treatment Effects (R-Learner)",
     xlab = "Estimated CATE",
     col = "lightblue")

#define function to plot marginal effects
plot_effect <- function(covar, label) {
  ggplot(data.frame(x = covar, y = cate_estimates), aes(x, y)) +
    geom_point(alpha = 0.3) +
    geom_smooth(method = "loess", color = "red") +
    labs(x = label, y = "Estimated CATE") +
    theme_minimal()
}

# Generate individual marginal plots
plot_effect(X$age, "Age")
plot_effect(X$inc, "Income")
plot_effect(X$educ, "Education")

# Estimate average treatment effect and confidence interval
avg_cate <- mean(cate_estimates, na.rm = TRUE)
se_cate <- sd(cate_estimates, na.rm = TRUE) / sqrt(sum(!is.na(cate_estimates)))
cat(sprintf("Average CATE: %.2f\n95%% CI: [%.2f, %.2f]\n",
            avg_cate, avg_cate - 1.96 * se_cate, avg_cate + 1.96 * se_cate))

# check for missing predictions (should be zero)
mean(is.na(cate_estimates))  

# Estimate and assess propensity scores
ps_model <- SuperLearner(Y = A, X = X, family = binomial(),
                         SL.library = c("SL.mean", "SL.glm"))  # Same library for consistency
ps <- predict(ps_model, X, onlySL = TRUE)$pred

# display summary of propensity scores
summary(ps)

# Plot distribution of estimated propensity scores
hist(ps,
     main = "Distribution of Propensity Scores",
     xlab = "Estimated P(Treatment | X)",
     col = "lightgreen", breaks = 30)

#create and label combined plot layout
p1 <- ggplot(data.frame(x = cate_estimates), aes(x)) +
  geom_histogram(fill = "lightblue", bins = 30) +
  labs(x = "CATE", y = "Frequency") +
  annotate("text", x = Inf, y = Inf, label = "(a)", hjust = 1.2, vjust = 1.2, size = 5)

p2 <- ggplot(data.frame(x = X$age, y = cate_estimates), aes(x, y)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", color = "red") +
  labs(x = "Age", y = "CATE") +
  annotate("text", x = Inf, y = Inf, label = "(b)", hjust = 1.2, vjust = 1.2, size = 5)

p3 <- ggplot(data.frame(x = X$inc, y = cate_estimates), aes(x, y)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", color = "red") +
  labs(x = "Income", y = "CATE") +
  annotate("text", x = Inf, y = Inf, label = "(c)", hjust = 1.2, vjust = 1.2, size = 5)

p4 <- ggplot(data.frame(x = X$educ, y = cate_estimates), aes(x, y)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", color = "red") +
  labs(x = "Education", y = "CATE") +
  annotate("text", x = Inf, y = Inf, label = "(d)", hjust = 1.2, vjust = 1.2, size = 5)

# display all plots in a 2x2 grid layout
(p1 | p2) / (p3 | p4)


############### Task 5 #######################################################################################################
# load and prepare data
data(pension)
Y <- pension$net_tfa
A <- pension$e401
L_full <- pension[, c("age", "inc", "educ", "fsize", "marr", "twoearn", "db", "pira", "hown")]
L_focus <- pension[, c("age", "inc", "educ")]

# Cap outcome to ±$1M
Y <- pmin(pmax(Y, -1e6), 1e6)

#cross-fitting setup
set.seed(123)
K <- 5
folds <- sample(rep(1:K, length.out = nrow(pension)))

#storage
pseudo_Y1 <- numeric(nrow(pension))
Y1_hat <- numeric(nrow(pension))
all_test_ps <- numeric(nrow(pension))
all_weights <- numeric(nrow(pension))

#Mainloop
for (k in 1:K) {
  train_idx <- which(folds != k)
  test_idx <- which(folds == k)
  
  # Propensity scores using full covariates
  ps_model <- SuperLearner(Y = A[train_idx],
                           X = L_full[train_idx, ],
                           family = binomial(),
                           SL.library = c("SL.glm", "SL.ranger", "SL.glmnet"))
  test_ps <- pmax(pmin(predict(ps_model, newdata = L_full[test_idx, ], onlySL = TRUE)$pred, 0.9), 0.1)
  
  # Outcome model using only age, inc, educ
  treated_train_idx <- train_idx[A[train_idx] == 1]
  Y1_model <- SuperLearner(Y = Y[treated_train_idx],
                           X = L_focus[treated_train_idx, ],
                           SL.library = c("SL.glm", "SL.ranger", "SL.glmnet"))
  test_Y1 <- predict(Y1_model, newdata = L_focus[test_idx, ], onlySL = TRUE)$pred
  
  # pseudo-outcome
  weights <- A[test_idx] / test_ps
  test_pseudo <- test_Y1 + weights * (Y[test_idx] - test_Y1)
  test_pseudo <- pmin(pmax(test_pseudo, -1e6), 1e6)
  
  #Store
  pseudo_Y1[test_idx] <- test_pseudo
  Y1_hat[test_idx] <- test_Y1
  all_test_ps[test_idx] <- test_ps
  all_weights[test_idx] <- weights
}

# Prediction error
residuals <- pseudo_Y1 - Y1_hat
mse <- mean(residuals^2)
mad_error <- mad(residuals, constant = 1)

set.seed(123)
boot_ci <- replicate(500, mean(sample(residuals^2, replace = TRUE)))
ci_low <- quantile(boot_ci, 0.025)
ci_high <- quantile(boot_ci, 0.975)

#  diagonistics and output
hist(all_test_ps, breaks = 20, main = "Propensity Scores (Trimmed 0.10–0.90)",
     xlab = "P(A=1|L)", col = "lightblue")

hist(all_weights, breaks = 50, main = "Weight Distribution (A/ps)",
     xlab = "Weight Value", col = "lightgreen")

hist(pseudo_Y1, breaks = 50, main = "Predicted Y¹ (DR-learner)",
     xlab = "Dollars ($)", col = "lightpink")

cat("\n-- Quick Check --\n")
cat("PS range:", round(min(all_test_ps), 3), "-", round(max(all_test_ps), 3), "\n")
cat("Weights: min =", round(min(all_weights), 2), ", med =", round(median(all_weights), 2), ", max =", round(max(all_weights), 2), "\n")
cat("Weights > 10:", sum(all_weights > 10), "\n")
cat("Capped at 1M:", round(mean(pseudo_Y1 == 1e6) * 100, 2), "% | at -1M:", round(mean(pseudo_Y1 == -1e6) * 100, 2), "%\n")

cat("\n-- Results --\n")
cat("Y1: from $", round(min(pseudo_Y1)/1000, 1), "k to $", round(max(pseudo_Y1)/1000, 1), "k\n")
cat("MAD error: $", round(mad_error), "\n")
cat("MSE: $", round(mse / 1e6, 2), "M\n")
cat("95% CI: [", round(ci_low / 1e6, 2), "-", round(ci_high / 1e6, 2), "] M\n")

# combine 3 diagnostic plots into one figure

# Plot a) Predicted Y¹
p1 <- ggplot(data.frame(Y1 = pseudo_Y1), aes(x = Y1)) +
  geom_histogram(fill = "pink", color = "black", bins = 30) +
  labs(
    title = "a",
    x = "Predicted Y¹ ($)",
    y = "Frequency"
  ) +
  theme_minimal()

# Plot b) Propensity Scores
p2 <- ggplot(data.frame(PS = all_test_ps), aes(x = PS)) +
  geom_histogram(fill = "skyblue", color = "black", bins = 20) +
  labs(
    title = "b",
    x = "P(A=1|X)",
    y = "Frequency"
  ) +
  theme_minimal()

# Plot c) Weights
p3 <- ggplot(data.frame(W = all_weights), aes(x = W)) +
  geom_histogram(fill = "lightgreen", color = "black", bins = 30) +
  labs(
    title = "c",
    x = "Weight",
    y = "Frequency"
  ) +
  theme_minimal()

# combine them vertically
combined_plot <- p1 / p2 / p3

# Show the plot
print(combined_plot)

