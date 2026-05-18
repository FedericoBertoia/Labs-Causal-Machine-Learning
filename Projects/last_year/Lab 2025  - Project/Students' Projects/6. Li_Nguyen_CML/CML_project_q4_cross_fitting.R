# load data and package
library(hdm)
library(ggplot2)
library(SuperLearner)  # For ML algorithms
library(caret)
library(rlearner)

data(pension)
# Crossfitting note: normally, we should use 5-10 folds, but bc it takes alot of time to run on my computer => so I only use 2,3

## SuperLearner library
sl_lib <- c("SL.glm", "SL.ranger")
cv_control <- SuperLearner.CV.control(V = 5) # cross-validation

# # Define variables
Y <- pension$net_tfa  # outcome, net total financial assets
A <- pension$e401     # treatsment, eligibility for 401(k), binary (1 = eligible)
X <- pension[, c("age","inc","fsize","educ", # Covariates (confounders)
                 "marr","twoearn","db", "ira","hown")]

pseudo_all <- matrix(NA,nrow(pension),2)
ID_pseudo <- 1:nrow(pension)
pseudo_all <- cbind(pseudo_all,ID_pseudo)

##### # 5-fold sample splitting
# Sample splitting
set.seed(1234)
folds <- createFolds(A,k=5)

# Initialize storage
n <- length(Y)
pseudo_outcome <- rep(NA, n)
propensity_score <- rep(NA, n)
weights <- rep(NA, n)

# cross fitting on propensity score and Y_hat
for(fold in 1:5){
  
  idx_test <- folds[[fold]]
  idx_train <- setdiff(1:n, idx_test)
  
  # Split data
  X_train <- X[idx_train,]
  X_test <- X[idx_test,]
  A_train <- A[idx_train]
  Y_train <- Y[idx_train]
  A_test <- A[idx_test]
  Y_test <- Y[idx_test]
  
  # Estimate Nuisance Parameters
  # predict outcome E(Y|L)
  model_Y <- SuperLearner(Y = Y_train, X = X_train, SL.library = sl_lib)
  Y_hat <- predict(model_Y, newdata = X_test)$pred
  
  # predict propensity score
  model_p <- SuperLearner(Y = A_train, X = X_train, SL.library = sl_lib, family = binomial())
  p_hat <- predict(model_p, newdata = X_test)$pred
  
  # compute pseudo outcome
  pso <- (Y_test - Y_hat) / (A_test - p_hat) 
  
  # compute weight
  weight <- (A_test - p_hat)^2
  
  # Store results
  pseudo_outcome[idx_test] <- pso
  propensity_score[idx_test] <- p_hat
  weights[idx_test] <- weight
}

# Regress the pseudooutcomes using weighted superleaner
Y_tilde <- rep(NA, n)
# cross fitting on pseudo outcome
for(fold in 1:5){
  idx_test <- folds[[fold]]
  idx_train <- setdiff(1:n, idx_test)
  #Predict Counterfactual Net Assets 
  cf_model <- SuperLearner(Y = pseudo_outcome[idx_train], X = X[idx_train,c('age','educ','inc')], 
                           SL.library = sl_lib, obsWeights=weights[idx_train])
  Y_t <- predict(cf_model, newdata = X[idx_test,c('age','educ','inc')], onlySL = TRUE)$pred
  Y_tilde[idx_test] <- Y_t
}

# Add CATEs to data 
pension$CATE <- Y_tilde
# Estimate Prediction Error and 95% Confidence Interval
theta_hat <- mean(Y_tilde)
eif <- Y_tilde - theta_hat
# Standard error and 95% CI
se <- sd(eif) / sqrt(length(Y_tilde))
ci_low <- theta_hat - 1.96 * se
ci_high <- theta_hat + 1.96 * se
cat("Estimated prediction error:", theta_hat, "\n")
cat("95% CI: [", ci_low, ",", ci_high, "]\n")


# Histogram of ITEs age
hist(pension$CATE, xlab = "Estimated CATE (age, education, income", 
     main = "Distribution of Estimated Effects", 
     col = "skyblue", breaks = 40)

# AGE
Y_tilde_age <- rep(NA, n)
# cross fitting on pseudo outcome
for(fold in 1:5){
  idx_test <- folds[[fold]]
  idx_train <- setdiff(1:n, idx_test)
  #Predict Counterfactual Net Assets 
  cf_model <- SuperLearner(Y = pseudo_outcome[idx_train], X = X[idx_train,c('age'), drop=FALSE], 
                           SL.library = sl_lib, obsWeights=weights[idx_train])
  Y_t <- predict(cf_model, newdata = X[idx_test,c('age'), drop=FALSE], onlySL = TRUE)$pred
  Y_tilde_age[idx_test] <- Y_t
}


# INCOME
Y_tilde_income <- rep(NA, n)
# cross fitting on pseudo outcome
for(fold in 1:5){
  idx_test <- folds[[fold]]
  idx_train <- setdiff(1:n, idx_test)
  #Predict Counterfactual Net Assets 
  cf_model <- SuperLearner(Y = pseudo_outcome[idx_train], X = X[idx_train,c('inc'), drop=FALSE], 
                           SL.library = sl_lib, obsWeights=weights[idx_train])
  Y_t <- predict(cf_model, newdata = X[idx_test,c('inc'), drop=FALSE], onlySL = TRUE)$pred
  Y_tilde_income[idx_test] <- Y_t
}


# EDUCATION
Y_tilde_edu <- rep(NA, n)
# cross fitting on pseudo outcome
for(fold in 1:5){
  idx_test <- folds[[fold]]
  idx_train <- setdiff(1:n, idx_test)
  #Predict Counterfactual Net Assets 
  cf_model <- SuperLearner(Y = pseudo_outcome[idx_train], X = X[idx_train,c('educ'), drop=FALSE], 
                           SL.library = sl_lib, obsWeights=weights[idx_train])
  Y_t <- predict(cf_model, newdata = X[idx_test,c('educ'), drop=FALSE], onlySL = TRUE)$pred
  Y_tilde_edu[idx_test] <- Y_t
}


# AGE
# Add CATEs to data 
pension$CATE_age <- Y_tilde_age
# Estimate Prediction Error and 95% Confidence Interval
theta_hat_age <- mean(Y_tilde_age)
eif_age <- Y_tilde_age - theta_hat_age
# Standard error and 95% CI
se_age <- sd(eif_age) / sqrt(length(Y_tilde_age))
ci_low_age <- theta_hat_age - 1.96 * se_age
ci_high_age <- theta_hat_age + 1.96 * se_age
cat("Estimated prediction error (Age):", theta_hat_age, "\n")
cat("95% CI (Age): [", ci_low_age, ",", ci_high_age, "]\n")

# Histogram of ITEs age
hist(pension$CATE_age, xlab = "Estimated CATE (age)", 
     main = "Distribution of Estimated Effects", 
     col = "skyblue", breaks = 40)
# scatter plot
ggplot(pension, aes(x = age, y = CATE_age)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", color = "blue") +
  theme_minimal() +
  labs(title = "CATE vs Age", x = "Age", y = "Estimated Effect")

# INCOME
# Add CATEs to data
pension$CATE_inc <- Y_tilde_income
# Estimate Prediction Error and 95% Confidence Interval
theta_hat_inc <- mean(Y_tilde_income)
eif_inc <- Y_tilde_income - theta_hat_inc
# Standard error and 95% CI
se_inc <- sd(eif_inc) / sqrt(length(Y_tilde_income))
ci_low_inc <- theta_hat_inc - 1.96 * se_inc
ci_high_inc <- theta_hat_inc + 1.96 * se_inc
cat("Estimated prediction error (Income):", theta_hat_inc, "\n")
cat("95% CI (Income): [", ci_low_inc, ",", ci_high_inc, "]\n")
# Histogram of ITEs income
hist(pension$CATE_inc, xlab = "Estimated CATE (income)", 
     main = "Distribution of Estimated Effects", 
     col = "skyblue", breaks = 40)
# by income
ggplot(pension, aes(x = inc, y = CATE_inc)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", color = "green") +
  theme_minimal() +
  labs(title = "CATE vs Income", x = "Income", y = "Estimated Effect")

# EDUCATION
# Add CATEs to data
pension$CATE_educ <- Y_tilde_edu
# Estimate Prediction Error and 95% Confidence Interval
theta_hat_edu <- mean(Y_tilde_edu)
eif_edu <- Y_tilde_edu - theta_hat_edu
# Standard error and 95% CI
se_edu <- sd(eif_edu) / sqrt(length(Y_tilde_edu))
ci_low_edu <- theta_hat_edu - 1.96 * se_edu
ci_high_edu <- theta_hat_edu + 1.96 * se_edu
cat("Estimated prediction error (Education):", theta_hat_edu, "\n")
cat("95% CI (Education): [", ci_low_edu, ",", ci_high_edu, "]\n")
# Histogram of ITEs income
hist(pension$CATE_educ, xlab = "Estimated CATE (education)", 
     main = "Distribution of Estimated Effects", 
     col = "skyblue", breaks = 40)
# by education
ggplot(pension, aes(x = educ, y = CATE_educ)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", color = "purple") +
  theme_minimal() +
  labs(title = "CATE vs Years of Education", 
       x = "Education", y = "Estimated Effect")



