# load data and package
library(hdm)
library(ggplot2)
library(SuperLearner)  # For ML algorithms
library(caret)

data(pension)

## SuperLearner library
sl_lib <- c("SL.glm", "SL.ranger", "SL.glmnet")
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
model_Y1s <- rep(NA,n)

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
  # Fit E[Y|A=1,X] and E[Y|A=0,X]
  model_Y1 <- SuperLearner(Y = Y_train[A_train == 1], X = X_train[A_train == 1,], 
                           SL.library=sl_lib, cvControl = cv_control)
  
  Y1_hat <- predict(model_Y1, newdata = X_test, OnlySL = TRUE)$pred # prediction come from ensemble
  
  # Propensity Score:
  # Estimate P(A|X):

  prop_model <- SuperLearner(Y = A_train, X = X_train, family = binomial(),
                             SL.library = sl_lib, cvControl = cv_control)
  p_hat <- predict(prop_model, newdata = X_test, OnlySL = TRUE)$pred
  
  # Compute Pseudo-Outcomes for DR-Learner
  Ci <- Y1_hat + (A_test / p_hat) * (Y_test - Y1_hat)
  
  # Store results
  pseudo_outcome[idx_test] <- Ci
  propensity_score[idx_test] <- p_hat
  model_Y1s[idx_test] <- model_Y1
}

Y1_tilde <- rep(NA, n)

# cross fitting on pseudo outcome
for(fold in 1:5){
  idx_test <- folds[[fold]]
  idx_train <- setdiff(1:n, idx_test)
  #Predict Counterfactual Net Assets 
  cf_model <- SuperLearner(Y = pseudo_outcome[idx_train], X = X[idx_train,c("age","inc","educ")], 
                           SL.library = sl_lib, cvControl = cv_control)
  Y1_t <- predict(cf_model, newdata = X[idx_test,c("age","inc","educ")], onlySL = TRUE)$pred
  Y1_tilde[idx_test] <- Y1_t
}

#plot histogram
hist(Y1_tilde, xlab="Predict counterfactual net assets (eligibility for 401k) ", main="Histogram of Counterfactual Net Assets (eligibility for 401k)")

model <- SuperLearner(Y = Y[A == 1], X = X[A == 1,c("age","inc","educ")], 
                                   SL.library=sl_lib, cvControl = cv_control)
Y1_hats <- predict(model, newdata = X[,c("age","inc","educ")], onlySL = TRUE)$pred

# Estimate Prediction Error and 95% Confidence Interval

Z <- (Y - Y1_tilde)^2

EZ_1s <- rep(NA, n)

# cross fitting on pseudo outcome
for(fold in 1:5){
  idx_test <- folds[[fold]]
  idx_train <- setdiff(1:n, idx_test)
  Z_train <- Z[idx_train]
  A_train <- A[idx_train]
  X_train <- X[idx_train,c("age","inc","educ")]
  #Predict E[Z|A=1,L]
  Z_model <- SuperLearner(Y = Z_train[A_train==1], X = X_train[A_train==1,], 
                           SL.library = sl_lib, cvControl = cv_control)
  EZ_1 <- predict(Z_model, newdata = X[idx_test,c("age","inc","educ")], onlySL = TRUE)$pred
  EZ_1s[idx_test] <- EZ_1
}

# compute E[Y^1] by ipw
EY1_ipw <- mean((A * EZ_1s) / propensity_score)

eif <- A/propensity_score*(Z-EZ_1s) + EZ_1s - EY1_ipw


# Standard error and 95% CI
se <- sd(eif) / sqrt(length(pension))
ci_low <- mean(eif) - 1.96 * se
ci_high <- mean(eif) + 1.96 * se

cat("Estimated prediction error:", mean(eif), "\n")

cat("95% CI: [", ci_low, ",", ci_high, "]\n")

#Estimated prediction error: 3040158871 
#95% CI: [ -15780598804 , 21860916546 ]

# checking for extreme influence
## 1. Check for Large Pseudo-Outcomes (Ci)
# Basic summary statistics
summary(pseudo_outcome)

# Histogram or boxplot
hist(pseudo_outcome, breaks = 50, main = "Pseudo-outcomes (Ci)", xlab = "Ci")
boxplot(pseudo_outcome, horizontal = TRUE, main = "Boxplot of Ci")


#There are many values that are far from the bulk of the data (e.g., hundreds of thousands away from the median), they may be influential points.


#2. Influence Function Approximation
plot(eif, main = "Residuals for Treated Units", ylab = "Residual", xlab = "Index")
abline(h = 0, col = "red")
# number of  potential extremes
length(which(abs(eif) > 3 * sd(eif)))
# 16

#3. Compare Predictions with Observed Outcomes - large discrepancies can be a signal
plot(Y1_tilde[A == 1], Y[A == 1], xlab = "Predicted Y1", ylab = "Observed Y", 
     main = "Observed vs Predicted (Eligible Units)")
abline(a = 0, b = 1, col = "red")

# 4. Check Influence from Propensity Score Weights
#Extremely small propensity scores (close to 0 or 1) cause large weights in pseudo-outcomes

summary(propensity_score) # Check for near-zero or near-one values
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.07297 0.21756 0.35120 0.37240 0.50392 0.91191 

hist(p_hat, main = "Histogram of Estimated Propensity Scores", xlab = "p_hat")
# Check instability of weights
weights <- A / propensity_score
boxplot(weights, main = "Boxplot of Inverse Propensity Weights")
#Extremely small propensity scores (close to 0 or 1) cause large weights in pseudo-outcomes

