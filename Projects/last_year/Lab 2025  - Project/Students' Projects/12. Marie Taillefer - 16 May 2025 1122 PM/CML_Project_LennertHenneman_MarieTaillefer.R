###################################################
################ CML Project 2025 #################
######## Lennert Henneman, Marie Taillefer ########
###################################################

### load packages ------------------------------
library(ggplot2)
library(SuperLearner)
library(npcausal)
library(tmle)
library(latex2exp)

### Load data ---------------------------------
library(hdm)
data(pension)
Y <- pension$net_tfa
A <- pension$e401
N <- length(A)
X <- pension[, c("age", "inc", "fsize", "educ", "marr", "twoearn", "db",
                 "pira", "hown")]
SL_lib <- c("SL.glm", "SL.glmnet", "SL.randomForest", "SL.earth")


##### Exercise 2 --------------------------------

##                Part 1: AIPW                 ##
ATE_AIPW <- ate(y = Y, a = A, x = X, nsplits = 5, sl.lib = SL_lib)

# Result examination
g1 <- ATE_AIPW$nuis$pi_a1
Q1 <- ATE_AIPW$nuis$mu_a1
Q0 <- ATE_AIPW$nuis$mu_a0
Q <-  Q1 * A +  Q0 * (1 - A)
EIC <- ((A / g1) - ( 1 - A) / (1 - g1))*(Y - Q) + Q1 - Q0 - 8091.754

min(g1) 
max(g1) 

AIPW_data <- data.frame(
  g1 = g1,
  A = factor(A),  # ensure A is treated as a categorical variable
  EIC = EIC
)

# examine propensity scores
ggplot(data = AIPW_data) + 
  geom_histogram(mapping = aes(x = g1, fill = A),
                 position = 'identity', alpha = 0.5)

# examine EICs
ggplot(data = AIPW_data) + 
  geom_boxplot(mapping = aes(x = EIC) )



##                Part 2: TMLE                 ##

Q_df <- cbind(Q0, Q1)
tmle(Y, A, X, Q = Q_df, g1W = g1)



##        Part 3: No covariate adjustment      ##

# manual computation

#  3. 1: ATE:
N_folds <- 5
set.seed(129)
inds <- sample(rep(seq_len(N_folds), length.out = N))

ATE_ests_nocov <- vector(length = N_folds)
IC_ests_nocov <- propensity <- vector(length = N)

Q_nocov_tot <- g1_nocov_tot <- Q1_nocov_tot <- Q0_nocov_tot <- vector(length = N)


for (fold in seq_len(N_folds)) {
  train <- inds != fold
  test <- inds == fold
  cat("computing ATE for fold: ", fold, "\n")
  
  # compute ATE for no adjustment
  
  propensity_nocov <- mean(A[train])
  g1_nocov <- propensity_nocov
  g0_nocov <- 1 - propensity_nocov
  Q1_nocov <- mean(Y[(A == 1) & train])
  Q0_nocov <- mean(Y[(A == 0) & train])
  Q_nocov <- ifelse(A[test] == 1, Q1_nocov, Q0_nocov)
  
  ate_aipw_nocov <- mean(((A[test] / g1_nocov) - (1 - A[test]) / g0_nocov) * (Y[test] - Q_nocov) + Q1_nocov - Q0_nocov)
  influence_curve_aipw_nocov <- ((A[test] / g1_nocov) - (1 - A[test]) / g0_nocov) * (Y[test] - Q_nocov) + Q1_nocov - Q0_nocov - ate_aipw_nocov
  
  ATE_ests_nocov[fold] <- ate_aipw_nocov
  IC_ests_nocov[test] <- influence_curve_aipw_nocov
  g1_nocov_tot[test] <- g1_nocov
  Q_nocov_tot[test] <- Q_nocov
  Q1_nocov_tot[test] <- Q1_nocov
  Q0_nocov_tot[test] <- Q0_nocov 
}

ate_aipw_nocov <- mean(ATE_ests_nocov)
ate_aipw_se_nocov <- 1/sqrt(N) * sd(IC_ests_nocov)

# print results
cat("ATE AIPW = ", ate_aipw_nocov)
cat("standard deviation = ", ate_aipw_se_nocov)
cat(" 95% CI = [", ate_aipw_nocov - 1.96 * ate_aipw_se_nocov,
    ", ", ate_aipw_nocov + 1.96 * ate_aipw_se_nocov, "]"  )


# 3.2: TMLE

# rescale Y (following (Frank HA, Karim ME, 2023))
min_y <- min(Y)
max_y <- max(Y)
scale <- max_y - min_y
Ys <- (Y - min_y) / scale

N_folds <- 5
set.seed(131)
inds <- sample(rep(seq_len(N_folds), length.out = N))

ATE_ests_nocov_tmle <- vector(length = N_folds)
IC_ests_nocov_tmle <- propensity <- vector(length = N)

for (fold in seq_len(N_folds)) {
  train <- inds != fold
  test <- inds == fold
  cat("computing ATE for fold: ", fold, "\n")
  
  # compute ATE for no adjustment
  
  propensity_nocov <- mean(A[train])
  g1_nocov <- rep(propensity_nocov, sum(test))
  g0_nocov <- 1 - g1_nocov
  Q1_nocov <- rep(mean(Ys[(A == 1) & train]), sum(test))
  Q0_nocov <- rep(mean(Ys[(A == 0) & train]), sum(test))
  Q_nocov <- ifelse(A[test] == 1, Q1_nocov, Q0_nocov)
  
  # clever covariate
  H1 <- A / g1_nocov
  H0 <- (1 - A)/ g0_nocov
  
  H <- H1-H0
  
  # Estimate parameters delta in the fluctuation models
  delta1_tmle <- coef(glm(Ys[test] ~ -1 + offset(qlogis(Q1_nocov)) + H1[test],
                          family = 'quasibinomial'))
  delta0_tmle <- coef(glm(Ys[test] ~ -1 + offset(qlogis(Q0_nocov)) + H0[test],
                          family = 'quasibinomial'))
  
  # Targeting step, scale back
  Q1_tmle <- plogis(qlogis(Q1_nocov) + delta1_tmle / g1_nocov) * scale + min_y
  Q0_tmle <- plogis(qlogis(Q0_nocov) + delta0_tmle / g0_nocov) * scale + min_y
  Q_tmle <- ifelse(A[test] == 1, Q1_tmle, Q0_tmle) 
  
  # Compute E[Y^1] and E[Y^0]
  EY1_tmle <- mean(Q1_tmle)
  EY0_tmle <- mean(Q0_tmle)
  ATE_tmle <- EY1_tmle - EY0_tmle
  
  influence_curve_tmle_nocov <- ((A[test] / g1_nocov) - (1 - A[test]) / g0_nocov) * 
    (Y[test] - Q_tmle) + Q1_tmle - Q0_tmle - ATE_tmle
  
  
  ATE_ests_nocov_tmle[fold] <- ATE_tmle
  IC_ests_nocov_tmle[test] <- influence_curve_tmle_nocov
  
}

ate_tmle_nocov <- mean(ATE_ests_nocov_tmle)
ate_tmle_se_nocov <- 1/sqrt(N) * sd(IC_ests_nocov_tmle)

cat("ATE = ", ate_tmle_nocov)
cat("standard deviation = ", ate_tmle_se_nocov)
cat(" 95% CI = [", ate_tmle_nocov - 1.96 * ate_tmle_se_nocov,
    ", ", ate_tmle_nocov + 1.96 * ate_tmle_se_nocov, "]"  )





# Exercise 3 ----------------------------------


##      Part 1: Calculate Pseudo Y       ##


# a function that calculates theta for each cross validation split.
Compute_theta_fold <- function(Y, A, X, train, test, sl_lib){
  
  # step 1: compute propensity estimates
  model <- SuperLearner(A[train], X[train, ], newX = X[test, ],
                        family=binomial(link = "logit"),
                        SL.library = sl_lib)
  g1 <- model$SL.predict # propensity estimates
  cat("propensity score calculation complete \n")
  print(model$coef)
  cat("\n")
  
  # step 2: compute counterfactual estimates
  model0 <- SuperLearner(Y[(A == 0) & train], X[(A == 0) & train, ],
                         newX = X[test, ], family = gaussian(), SL.library = sl_lib)
  Q0 <- model0$SL.predict
  cat("Q0 calculation complete \n")
  print(model0$coef)
  cat("\n")
  
  model1 <- SuperLearner(Y[(A == 1) & train], X[(A == 1) & train, ],
                         newX = X[test, ], family = gaussian(), SL.library = sl_lib)
  Q1 <- model1$SL.predict
  cat("Q1 calculation complete \n")
  print(model1$coef)
  cat("\n")
  Q <- A[test] * Q1 + (1 - A[test]) * Q0
  
  # step 3: calculate Theta
  theta <- mean((A[test] - g1) * (Y[test] - Q))
  influence_curve_theta <- (A[test] - g1) * (Y[test] - Q) - theta
  
  return(list(theta = theta, IC = influence_curve_theta))
}



# Now splitting the data into folds again and computing theta
N <- length(A)
N_folds <- 5
set.seed(131)
inds <- sample(rep(seq_len(N_folds), length.out = N))

Theta_ests <- vector(length = N_folds)
IC_ests_theta <- propensity <- vector(length = N)


for (fold in seq_len(N_folds)) {
  train <- inds != fold
  test <- inds == fold
  cat("computing Theta for fold: ", fold, "\n")
  
  output <- Compute_theta_fold(Y, A, X, train, test, SL_lib)
  
  Theta_ests[fold] <- output$theta
  IC_ests_theta[test] <- output$IC
  
}

Theta <- mean(Theta_ests)
Theta_sd <- 1/sqrt(N) * sd(IC_ests_theta)
# p-value:
z <- (0 - Theta) / Theta_sd 
p <- (1 - pnorm(abs(z))) * 2 

cat("Theta = ", Theta)
cat("standard deviation = ", Theta_sd)
cat(" 95% CI = [", Theta - 1.96 * Theta_sd,
    ", ", Theta + 1.96 * Theta_sd, "]"  )
cat("Z-score = ", z)
cat("p-value =", p)

# Exercise 4 ------------------------------------


detach("package:npcausal", unload = TRUE, character.only = TRUE)
library(SuperLearner)
library(pROC)
#Change SL.randomForest to SL.ranger for faster results
SL_lib <- c("SL.glm", "SL.glmnet", "SL.ranger", "SL.earth")

#Create a function for computing nuisance parameters using cross-fitting
cf_nuisance_param <- function(data, SL_lib, N_folds, seed = 123){
  
  set.seed(seed)
  
  Y <- data$net_tfa
  A <- data$e401
  X <- data[, c("age", "inc", "fsize", "educ", "marr",
                "twoearn", "db", "pira", "hown")]
  
  
  #Split data into folds
  
  N <- nrow(data)
  
  inds <- sample(rep(1:N_folds, length.out = N))
  
  data$Y_pred <- NA
  data$prop_score <- NA
  
  #Implement cross-fitting for nuisance parameters
  for (fold in 1:N_folds) {
    train <- inds != fold
    test <- inds == fold
    
    cat("\n", "Computing outcome predictions and propensity scores with cross-fitting.", "\n") 
    cat("Fold: ", fold, "\n")
    
    # Compute outcome predictions
    model_out_preds <-  SuperLearner(Y = Y[train], X = X[train,], SL.library = SL_lib)
    data$Y_pred[test] <- as.vector(predict(model_out_preds, 
                                           newdata = X[test,], 
                                           onlySL = TRUE)$pred)
    
    cat("\n", "Learner coefficients for outcome prediction model.", "\n")
    print(model_out_preds$coef)
    cat("Summary of outcome predictions", "\n")
    print(summary(data$Y_pred[test]))
    
    #Compute propensity scores
    model_prop_scores <- SuperLearner(Y = A[train], X = X[train,], 
                                      family = binomial(), 
                                      SL.library = SL_lib)
    
    data$prop_score[test] <- as.vector(predict(model_prop_scores, 
                                               newdata = X[test,], 
                                               onlySL = TRUE)$pred)
    
    cat("\n", "Learner coefficients for propensity score model.", "\n")
    print(model_prop_scores$coef)
    cat("Summary of propensity scores.", "\n")
    print(summary(data$prop_score[test]))
    
  }
  data
}

#Implement cross-fitting loop for rlearner

set.seed(123)
#Split data into five folds

N_folds <- 5

pension$rlearner_pred = NA #Initialiaze a variable to save predictions

inds <- sample(rep(1:N_folds, length.out = N))


for (fold in 1:N_folds) {
  train <- inds != fold
  test <- inds == fold
  
  cat("\n", "computing CATE for fold: ", fold, "\n")
  
  #Compute outcome predictions (based on L)
  pension_nuisance <- cf_nuisance_param(data = pension[train,], 
                                        SL_lib = SL_lib, N_folds = 5)
  
  
  #Quality check of outcome predictions for current fold
  cat( "\n", "Outcome fit (fold", fold, "):\n")
  mse_m  <- mean((Y[train] - pension_nuisance$Y_pred)^2)
  r2_m   <- 1 - mse_m / var(Y[train])
  cat("  Outcome fit: CV‑MSE =", round(mse_m,2),
      "| CV‑R2 =", round(r2_m,2), "\n")
  
  
  #Quality check of propensity scores for current fold
  par(mar = c(4, 4, 2, 1))  
  
  cat("\n", "Propensity scores (fold", fold, "):\n")
  print(summary(pension_nuisance$prop_score))
  
  plot_df <- data.frame(
    prop_score = pension_nuisance$prop_score,
    treatment = factor(A[train])
  )
  
  print(ggplot(plot_df, aes(x = prop_score, fill = treatment)) +
    geom_histogram(alpha = 0.5, position = "identity", bins = 30) +
    labs(title = paste("Fold", fold, "Propensity Scores by Treatment"),
         x = "Propensity Score", fill = "Treatment") +
    xlim(0, 1) +
    theme_minimal())
  
  
  n_low  <- sum(pension_nuisance$prop_score < 0.05)
  n_high <- sum(pension_nuisance$prop_score > 0.95)
  n_obs  <- nrow(pension_nuisance)
  
  
  cat("\n", "Extreme propensity scores (fold", fold, "):\n")
  cat(
    sprintf(
      "  Low (<0.05): %d (%.1f%%) | High (>0.95): %d (%.1f%%)\n",
      n_low,  100 * n_low / n_obs,
      n_high, 100 * n_high / n_obs
    )
  )
  
  auc_value <- auc(roc(response = A[train], predictor = pension_nuisance$prop_score))
  cat(sprintf("AUC for propensity model (fold %d): %.3f\n", fold, auc_value))
  
  
  # Compute the pseudooutcomes
  pension_nuisance$pseudooutcome <- 
    (Y[train] - pension_nuisance$Y_pred) / 
    (A[train] - pension_nuisance$prop_score)
  
  
  print(ggplot(data = pension_nuisance) + 
    geom_boxplot(mapping = aes(x = pseudooutcome)))
  
  
  # Regress the pseudooutcomes using weighted SuperLearner
  pension_nuisance$weights <- (A[train] - pension_nuisance$prop_score)^2
  model_Rlearner <- SuperLearner(Y = pension_nuisance$pseudooutcome, 
                                 X = X[train, c("age", "inc", "educ")], 
                                 SL.library = SL_lib, 
                                 obsWeights = pension_nuisance$weights)
  
  
  # Predict CATE on current fold
  pension$rlearner_pred[test] <- predict(model_Rlearner, 
                                         newdata = X[test, c("age", "inc", "educ")], 
                                         onlySL = TRUE)$pred
}


# Calculate ATE and 95%CI for sanity check

tau_hat <- pension$rlearner_pred

CATE_hat <- mean(tau_hat)

se_hat  <- sd(tau_hat) / sqrt(N)

ci_low  <- CATE_hat - 1.96 * se_hat
ci_high <- CATE_hat + 1.96 * se_hat

cat("ATE estimate =", round(CATE_hat,2), "\n")
cat("95% CI      = [", round(ci_low,2), ", ", round(ci_high,2), "]\n")

hist(pension$rlearner_pred, 
     main = "Distribution of Estimated Conditional Average Treatment Effects", 
     xlab = "Estimated CATE", 
     col = "skyblue", 
     border = "white", 
     breaks = 100)

ggplot(subset(pension, abs(rlearner_pred) <= 1e6), aes(x = age, y = rlearner_pred)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess", se = TRUE, col = "blue") +
  labs(title = "Treatment Effect by Age",
       x = "Age", y = "Estimated CATE")

ggplot(subset(pension, abs(rlearner_pred) <= 1e6), aes(x = inc, y = rlearner_pred)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess", se = TRUE, col = "darkgreen") +
  labs(title = "Treatment Effect by Income",
       x = "Income", y = "Estimated CATE")

ggplot(subset(pension, abs(rlearner_pred) <= 1e6), aes(x = educ, y = rlearner_pred)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess", se = TRUE, col = "darkred") +
  labs(title = "Treatment Effect by Education",
       x = "Years of Education", y = "Estimated CATE")




# Exercise 5 -------------------------------------

detach("package:pROC", unload = TRUE, character.only = TRUE)
SL_lib <- c("SL.glm", "SL.glmnet", "SL.randomForest", "SL.earth")


# setup function to compute pseudo Y for each fold
Compute_pseudo_y1_fold <- function(Y, A, X, train, test, sl_lib){
  
  # step 1: divide the training data in two groups:
  inds <- sample(rep(seq_len(2), length.out = sum(train)))
  Y_trainA  <- Y[train][inds == 1]
  Y_trainB  <- Y[train][inds == 2]
  A_trainA  <- A[train][inds == 1]
  A_trainB  <- A[train][inds == 2]
  Y1_trainA <- Y_trainA[A_trainA == 1]
  Y1_trainB <- Y_trainB[A_trainB == 1]
  X_trainA  <- X[train, ][inds == 1, ]
  X_trainB  <- X[train, ][inds == 2, ]
  X1_trainA <- X_trainA[A_trainA == 1, ]
  X1_trainB <- X_trainB[A_trainB == 1, ]
  
  
  # step 2: compute propensity estimates
  modelA <- SuperLearner(A_trainA, X_trainA, newX = X_trainB,
                         family=binomial(link = "logit"),
                         SL.library = sl_lib)
  
  cat("propensity trained on first half of train data \n")
  print(modelA$coef)
  cat("\n")
  
  modelB <- SuperLearner(A_trainB, X_trainB, newX = X_trainA,
                         family=binomial(link = "logit"),
                         SL.library = sl_lib)
  
  g1 <- vector(length = sum(train))
  g1[inds == 1] <- modelB$SL.predict
  g1[inds == 2] <- modelA$SL.predict
  
  cat("propensity score calculation complete \n")
  print(modelB$coef)
  cat("\n")
  
  
  # step 3: compute counterfactual estimates
  model1A <- SuperLearner(Y1_trainA, X1_trainA, newX = X_trainB,
                          family = gaussian(), SL.library = sl_lib)
  
  cat("Q1 trained on first half of train data \n")
  print(model1A$coef)
  cat("\n")
  
  model1B <- SuperLearner(Y1_trainB, X1_trainB, newX = X_trainA,
                          family = gaussian(), SL.library = sl_lib)
  
  Q1 <- vector(length = sum(train))
  Q1[inds == 1] <- model1B$SL.predict
  Q1[inds == 2] <- model1A$SL.predict
  
  cat("Q1 calculation complete \n")
  print(model1B$coef)
  cat("\n")
  
  
  # step 4: Target individual pseudo Y1 and predict it
  pseudo_Y1 <- Q1 + A[train] / g1 * (Y[train] - Q1)
  model2 <- SuperLearner(pseudo_Y1, X[train, ],
                         newX = X[test, ], family = gaussian(), SL.library = sl_lib)
  C <- model2$SL.predict
  cat("C calculation complete \n")
  print(model2$coef)
  cat("\n")
  
  
  return(list(C = C, propensity = g1, pseudo_Y1 = pseudo_Y1, Q1 = Q1))
}


# setup the cross fitting iterations
N_folds <- 5
set.seed(133)
inds <- sample(rep(seq_len(N_folds), length.out = N))

pseudo_Y1 <- prop2 <- DR_C <- vector(length = N)

for (fold in seq_len(N_folds)) {
  train <- inds != fold
  test <- inds == fold
  cat("computing C for fold: ", fold, "\n")
  
  output <- Compute_pseudo_y1_fold(Y, A, X, train, test, SL_lib)
  
  DR_C[test] <- output$C
  pseudo_Y1[test] <- output$pseudo_Y1
  prop2[test] <- output$propensity
  
}

# plot the results in a histogram
ggplot() + geom_histogram(aes(x = DR_C)) + 
  xlim(c(-0.5e5, 2e5)) +
  xlab(expression(tilde(Y)^1))



##            Part 2: estimate MSPE           ##


MSPE_AIPW <- ate(y = (Y - DR_C)^2, a = A, x = X, nsplits = 5, sl.lib = SL_lib)

# check predicted values
MSPE <- MSPE_AIPW$nuis$mu_a1
EIC <- A / g1 * ((Y - DR_C)^2 - MSPE) + (Y - DR_C)^2 - 3059733884
ggplot() + geom_boxplot(aes(x = MSPE))
ggplot() + geom_boxplot(aes(x = EIC))


