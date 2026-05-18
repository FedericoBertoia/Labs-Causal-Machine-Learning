#Garai Medrano Kareaga and Kylian Schalembier
#### LOADING LIBRARIES ###
#install.packages("earth")
#install.packages("randomForest")
#install.packages("ranger")
#install.packages("SuperLearner")
#install.packages("devtools")
#install_github("ehkennedy/npcausal")
#install_github("xnie/rlearner")
#install.packages("hdm")
#install.packages("tmle",dependencies = TRUE)
#install.packages("ICEbox")
library(ICEbox)
library(tmle)
library(ggplot2)
library(SuperLearner)
library(devtools)
library(npcausal)
library(hdm)
library(rlearner)
library(gbm)
data("pension")
head(pension)
attach(pension)
V <- 5
#################################
#################################
#########  EXERCISE 2  ##########
#################################
#################################

#### NAMING VARIABLES ###
Y <- pension$net_tfa #outcome
A <- pension$e401 #treatment          
X <- as.data.frame(pension[, c("age", "inc", "fsize", "educ", "marr", "pira", "twoearn", "db", "hown")]) #confounders

#### COMPUTING ATE/ATT W/ AIPW ###

sl_lib <- list("SL.glm","SL.step","SL.glm.interaction","SL.randomForest")
result <- ate(y = Y, a = A, x = X, nsplits = 5, sl.lib = sl_lib)
#parameter       est        se     ci.ll     ci.ul pval
#1      E{Y(0)} 14744.463 1133.0894 12523.607 16965.318    0
#2      E{Y(1)} 21921.426  853.9018 20247.779 23595.074    0
#3 E{Y(1)-Y(0)}  7176.963 1315.6724  4598.245  9755.681    0
result2 <- att(y = Y, a = A, x = X, nsplits = 5, sl.lib = sl_lib)
#parameter       est       se     ci.ll    ci.ul pval
#1      E(Y|A=1) 30347.389 1243.033 27911.045 32783.73    0
#2   E{Y(0)|A=1} 20689.203 1226.689 18284.892 23093.51    0
#3 E{Y-Y(0)|A=1}  9658.186 1523.856  6671.429 12644.94    0
#ATE is more suitable here because eligibility applies to the entire population, not just those who take it up. This informs policy-level decisions about offering eligibility to all.
influences<-result2$ifvals$V2-result2$ifvals$V3
influences_st<-influences/sd(influences)
max(influences)
max(influences_st)
sum(abs(influences)>3*sd(influences))
#Not too crazy

#### COMPUTING ATE/ATT W/ TMLE ###

sl_lib2 <- list("SL.glm","SL.step","SL.glm.interaction","SL.randomForest")

tmle_package <- tmle(Y = Y,
                     A = A,
                     W = X,
                     Q.SL.library = sl_lib2,
                     g.SL.library = sl_lib2, V.Q = 5, V.g = 5)
summary(tmle_package)

#Parameter Estimate:  7079.6
#Estimated Variance:  1759500
#p-value:  9.438e-08
#95% Conf Interval:  (4479.8, 9679.4)
#### NAIVE ESTIMATE ####
#mean(pension[pension$e401==1,]$net_tfa)-mean(pension[pension$e401==0,]$net_tfa)
naive=glm(Y~A)#duh
summary(naive)

#################################
#################################
#########  EXERCISE 3  ##########
#################################
#################################
# ----------------------------------------------------
# Test H0: theta = E{Cov(A, Y | L)} = 0
# ----------------------------------------------------
# Efficient influence curve:
#   phi_i(theta, eta) = (A_i - E[A | L_i]) * (Y_i - E[Y | L_i]) - theta
#   where eta = (E[A | L], E[Y | L]).
#
# Estimate E[A | L] and E[Y | L] with SuperLearner, use 5-fold cross-fitting.

print("Q3")

fold_id <- sample(rep(1:V, length.out = nrow(pension)))

theta_fold <- numeric(V)
infl <- numeric(nrow(pension))

for (v in 1:V) {
  train <- fold_id != v
  test <- !train
  
  # mu_hat(L) = E[Y | L]
  mu_fit <- SuperLearner(
    Y = Y[train],
    X = X[train, ],
    newX = X[test, ],
    SL.library = sl_lib
  )
  mu_hat <- mu_fit$SL.predict
  
  # pi_hat(L) = E[A | L]
  pi_fit <- SuperLearner(
    Y = A[train],
    X = X[train, ],
    newX = X[test, ],
    family = binomial(),
    SL.library = sl_lib
  )
  pi_hat <- pi_fit$SL.predict
  
  # fold-specific theta
  res_prod <- (A[test] - pi_hat) * (Y[test] - mu_hat)
  theta_fold[v] <- mean(res_prod)
  
  # influence values for variance
  infl[test] <- res_prod - theta_fold[v]
}

theta_hat <- mean(theta_fold)
se_hat <- sqrt(mean(infl^2) / nrow(pension))
z_stat <- theta_hat / se_hat
p_value <- 2 * (1 - pnorm(abs(z_stat)))

cat(
  "\n--- Question 3 results ---\n",
  "theta_hat =", round(theta_hat, 4), "\n",
  "SE        =", round(se_hat, 4), "\n",
  "z_stat    =", round(z_stat, 3), "\n",
  "p_value   =", signif(p_value, 3), "\n"
)

# Interpretation:
#   If p_value < 0.05, reject H0 -> evidence that eligibility
#   affects net assets via conditional covariance.
#   Otherwise, fail to reject H0.



#################################
#################################
#########  EXERCISE 4  ##########
#################################
#################################
#Yields a terrible model, honestly no idea why it doesn't work
X_control <- pension[,c("fsize", "marr", "twoearn", "db", "pira", "hown")]
X_hte <- as.data.frame(pension[,c("inc", "age", "educ")])  # Variables for HTE
sl_lib <- list("SL.ranger","SL.glm", "SL.gam", "SL.xgboost")
pension$e401<-as.factor(pension$e401)

#5-Fold Cross-Fitting Setup
folds <- sample(rep(1:5, length.out = nrow(pension)))
pension$cate <- NA  # To store CATEs
sl_lib <- c("SL.glm", "SL.step","SL.glm.interaction", "SL.mean")
#Main R-Learner Loop
fold=1
for (fold in 1:5) {
  train <- folds != fold
  test <- folds == fold
  train_df <- as.data.frame(cbind(Y = Y[train], A = A[train], X[train,]))
  test_df <- as.data.frame(cbind(A = A[test], X[test,]))
  if(sum(A[train] == 1) < 10 || sum(A[train] == 0) < 10) next
  # Step 1: Fit nuisance models
  Q_fit <- SuperLearner(
    Y = train_df$Y,
    X = train_df[, -1],
    SL.library = sl_lib,
    family = gaussian()
  )
  # Propensity model (g) - predicts treatment probability
  g_fit <- SuperLearner(
    Y = train_df$A,
    X = train_df[,3:ncol(train_df)],
    SL.library = sl_lib,
    family = binomial()
  )
  # Compute pseudo-outcomes
  Q_test <- predict(Q_fit, newdata = test_df)$pred
  g_test <- pmax(pmin(predict(g_fit, newdata = test_df[, -1])$pred, 0.99), 0.01)
  
  #R-learner equation
  pseudo_outcome <- (Y[test] - Q_test) / (A[test] - g_test) # TREATMENT IN DENOMINATOR
  pseudo_outcome <- pmin(pmax(pseudo_outcome, quantile(pseudo_outcome, 0.01, na.rm = TRUE)), 
                         quantile(pseudo_outcome, 0.99, na.rm = TRUE))
  weights <- (A[test] - g_test)^2 # SQUARED TREATMENT RESIDUALS
  valid_idx <- is.finite(pseudo_outcome) & weights > 1e-6
  if(sum(valid_idx) < 5) next 
  cate_data <- as.data.frame(X_hte[test,][valid_idx,])
  #Fit CATE model on effect modifiers
  cate_fit <- SuperLearner(
    Y = pseudo_outcome[valid_idx],
    X = cate_data, # Variables for HTE
    SL.library = list("SL.ranger","SL.xgboost","SL.gam","SL.earth", "SL.mean"),
    obsWeights = weights[valid_idx],
    family = gaussian()
  )
  cate_fit <- ranger(
    y = pseudo_outcome[valid_idx],
    x = cate_data,
    case.weights = weights[valid_idx],  # Use obsWeights
    num.trees = 500,
    mtry = max(floor(ncol(cate_data)/3)),importance = "none")
  summary(cate_fit)
  pension$cate[test][valid_idx] <- predict(cate_fit, cate_data)$pred
  print(fold)
}
# Alternative using rboost
X.rl=as.matrix(X)
A.rl=as.numeric(as.character(A))
Y.rl=as.numeric(Y)
rboost.fit <- rboost(X.rl , A.rl , Y.rl, k_folds=5)
summary(rboost.fit$p_hat)
summary(rboost.fit$weights)
rboost.fit$weights<-pmax(pmin(rboost.fit$weights, 0.99), 0.01)
rboost.est <- predict(rboost.fit, X.rl)
summary(rboost.est)
hist(rboost.est, 
     main = "Conditional Average Treatment Effect", 
     xlab = "Estimated treatment effect (USD)")
#Very weird

#Now, let's predict in function of income, age and education level.

cate_model_i <- ranger(rboost.est ~ inc, data = pension)
cate_model_a <- ranger(rboost.est ~ age, data = pension)
cate_model_e <- ranger(rboost.est ~ educ, data = pension)

# Create prediction data frames for each variable
pred_data_inc <- data.frame(inc = seq(min(pension$inc), max(pension$inc), length.out = 100))
pred_data_age <- data.frame(age = seq(min(pension$age), max(pension$age), length.out = 100))
pred_data_educ <- data.frame(educ = seq(min(pension$educ), max(pension$educ), length.out = 100))
# Get predictions from each model
pred_inc <- predict(cate_model_i, data = pred_data_inc)$predictions
pred_age <- predict(cate_model_a, data = pred_data_age)$predictions
pred_educ <- predict(cate_model_e, data = pred_data_educ)$predictions

# Create plots
par(mfrow = c(1, 3)) # Arrange plots side by side

# 1. Income plot
plot(pred_data_inc$inc, pred_inc, type = "l", 
     main = "Effect by income", 
     xlab = "Income", ylab = "Predicted effect")

# 2. Age plot
plot(pred_data_age$age, pred_age, type = "l", 
     main = "Effect by age", 
     xlab = "Age", ylab = "Predicted effect")

# 3. Education plot
plot(pred_data_educ$educ, pred_educ, type = "l", 
     main = "Effect by level of ducation", 
     xlab = "Education Level", ylab = "Predicted effect")

#################################
#################################
#########  EXERCISE 5  ##########
#################################
#################################


folds <- sample(1:V, nrow(pension), replace = TRUE)
Y1_hat <- numeric(nrow(pension))
k=1
for (k in 2:V) {
  train_idx <- folds != k
  test_idx <- folds == k
  
  # Train outcome model (A=1 only)
  Y_model <- SuperLearner(Y = Y[train_idx & A == 1], 
                          X = X[train_idx & A == 1, ], 
                          SL.library = sl_lib)
  Y1_train <- predict(Y_model, X[train_idx, ])$pred
  
  # Train propensity model
  p_model <- SuperLearner(Y = A[train_idx], 
                          X = X[train_idx, ], 
                          SL.library = sl_lib, 
                          family = binomial())
  p_train <- predict(p_model, X[train_idx, ])$pred
  
  # Compute pseudo-outcome C
  C_train <- Y1_train + (A[train_idx] / pmax(p_train, 0.05)) * (Y[train_idx] - Y1_train)
  
  # Train DR-learner
  dr_model <- SuperLearner(Y = C_train, 
                           X = X[train_idx, ], 
                           SL.library = sl_lib)
  Y1_hat[test_idx] <- predict(dr_model, X[test_idx, ])$pred
  print(k)
}
cl
hist(Y1_hat, breaks = 60, 
     main = "Distribution of Predicted Y¹",
     xlab = "Counterfactual Net Financial Assets",
     )
#Create the usual DR-learner model to get another model for Y^1
phat <- rboost.fit$p_hat
tboost.fit <- tboost(X.rl, A.rl, Y.rl)
y1.t <- tboost.fit$y_1_pred
y0.t <- tboost.fit$y_0_pred
ypseudo.drl <- (A.rl/phat)*(Y.rl-y1.t)+y1.t-(1-A.rl)/(1-phat)*(Y.rl-y0.t)-y0.t
drboost.fit<- gbm(ypseudo.drl~X.rl[,1]+X.rl[,2]+X.rl[,3]+X.rl[,4]+X.rl[,5]+X.rl[,6]+X.rl[,7]+X.rl[,8]+X.rl[,9],
                  distribution="gaussian", n.trees=5000, interaction.depth=4)
drboost.est<- predict(drboost.fit)

#Simulate the mean square error
boot_mse <- replicate(500, {
  idx <- sample(nrow(pension), replace = TRUE)
  (Y1_hat[idx] - y1.t[idx])^2
})
ci <- quantile(boot_mse, c(0.025, 0.975))
median(boot_mse)
ci
