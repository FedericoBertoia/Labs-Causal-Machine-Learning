"""
CML project R-script
"""

###set up
#Thomas
setwd("C:/Users/thoke/OneDrive/Documenten/1ste Master/semester 2/Causal Machine Learning/Project")
set.seed(123)

source("functions.R")

library(hdm) 
data(pension)
library(devtools)
library(npcausal)
library(SuperLearner)
library(earth)
library(ranger)
library(randomForest)
library(rlearner)
library(tmle)
?pension

#Florin
script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
# Set the working directory to the script's directory
setwd(script_dir)


set.seed(123)


# Loading data and other packages
install.packages("hdm")
library(hdm)
data(pension)

library(SuperLearner)
library(tidyverse)

###Q1



### QUESTION 2 (Thomas)
"""
Instructions:
Use AIPW as well as TMLE to estimate the effect of eligibility for enrolling in a 401(k) plan 
on net financial assets, controlling for age, income, family size, years of
education, a married indicator, a two-earner status indicator, a defined benefit pension status indicator, 
an IRA (individual retirement account) participation indicator, and a home ownership indicator. 
Report either the ATE or the ATT (or both), depending on what seems most suitable to you. 
Compare with the corresponding estimate when no covariate adjustment is made.
"""

A <- pension$e401      #eligibility (treatment)
Y <- pension$net_tfa   #net financial assets (outcome)
X <- pension[, c("age", "inc", "fsize", "educ", "marr", "twoearn", "db", "pira", "hown" )]   #covariates


sl_lib <- list("SL.glm","SL.step","SL.glm.interaction","SL.randomForest")
sl_lib3 <- list("glmnet", "SL.polymars", "SL.kernelKnn", "SL.glm")



#check propensity scores to decide on on att or ate
pihat_nonpar_fit <- estimate_pi_nonpar(A = A,
                                       X = X,
                                       sl_lib = sl_lib,
                                       nsplits = 5)

pihat_nonpar <- pihat_nonpar_fit$preds
pihat_nonpar_bounded <- bound_pi(pihat_nonpar)

# Check positivity assumption 
Apihat = as.data.frame(cbind(A, pihat_nonpar_bounded)) #here we are using estimated propensity scores with cross-fitting

pi_A1 <- Apihat$pihat_nonpar_bounded[Apihat$A == 1]
pi_A0 <- Apihat$pihat_nonpar_bounded[Apihat$A == 0]

pi_A1_inverse <- (1/Apihat$pihat_nonpar_bounded[Apihat$A == 1])
pi_A0_inverse <- (1/Apihat$pihat_nonpar_bounded[Apihat$A == 0])

plot(density(pi_A1, na.rm = TRUE), col = "red", lwd = 2, 
     main = "Kernel Density Estimation of pi_hat_nonpar_bounded", 
     xlab = "pi_hat_nonpar_bounded", ylim = c(0, max(density(pi_A1)$y, density(pi_A0)$y)))

lines(density(pi_A0, na.rm = TRUE), col = "blue", lwd = 2)

legend("topright", legend = c("A = 1", "A = 0"), col = c("red", "blue"), lwd = 2)

boxplot(pi_A1_inverse, pi_A0_inverse, 
        names = c("A = 1", "A = 0"), 
        col = c("red", "blue"), 
        main = "Boxplot of 1/pi_hat_nonpar_bounded", 
        ylab = "1/pi_hat_nonpar_bounded")




## AIPW part using np causal
#  att should be reported, as indicated by checking propensity scores

aipw_5cv <- att(y = Y, a = A, x = X, nsplits = 5, sl.lib = sl_lib2)
aipw_5cv
aipw_ate_5cv <- ate(y = Y, a = A, x = X, nsplits = 5, sl.lib = sl_lib)

#run again without covariates
sl_lib0 <- list("SL.mean")
X_empty <- pension[, c()]   #covariates
aipw_unadj <- att(y = Y, a = A, x = X_empty, nsplits = 5, sl.lib = sl_lib0)



##TMLE part

# without randomForest, but with SL.ranger instead 
# When including randomForest, R keeps aborting
sl_lib2 <- list("SL.glm","SL.step","SL.glm.interaction", "SL.ranger")

N        <- length(A)
nsplits  <- 5
split_id <- sample(rep(1:nsplits, ceiling(N/nsplits))[1:N])

Q0 <- rep(NA, N)
g <- rep(NA, N)

coefs.Q0 <- rep(0, length(sl_lib2))
coefs.g <- rep(0, length(sl_lib2))

for(v in 1:nsplits) {
  train <- split_id != v
  test  <- split_id == v
  
  Q0.sl <- SuperLearner(Y[(A == 0) & train], X[(A == 0) & train, ], newX = X[test, ], family = 'gaussian', SL.library = sl_lib2)
  Q0[test] <- Q0.sl$SL.predict
  coefs.Q0 <- coefs.Q0 + Q0.sl$coef
  
  g.sl <- SuperLearner(A[train], X[train, ], newX = X[test, ], family='binomial', SL.library = sl_lib2)
  g[test] <- g.sl$SL.predict
  coefs.g <- coefs.g + g.sl$coef
}


P_A1 <- mean(A)
g <- pmin(g, 0.99)
g <- pmax(g, 0.01)

# Construct clever covariate H
H <- ifelse(A == 1, 1 / P_A1, -g / (P_A1 * (1 - g)))

# Fit fluctuation model (update Q0 to Q0_star)
epsilon_model <- glm(Y ~ -1 + H + offset(Q0), family = gaussian())
epsilon <- coef(epsilon_model)["H"]
Q0_tmle <- Q0 + epsilon * H

# Compute ATT
theta_hat <- mean( (Y[A == 1] - Q0_tmle[A == 1]) )

# Influence function and standard error
IF_part1 <- (A / P_A1) * (Y - Q0_tmle - theta_hat)
IF_part2 <- ((1 - A) * g / (P_A1 * (1 - g))) * (Y - Q0_tmle)
IF <- IF_part1 - IF_part2
se <- sqrt(mean(IF^2)) / sqrt(N)

# Confidence interval
ci <- theta_hat + c(-1.96, 1.96) * se

# Results
results <- data.frame(
  ATT = theta_hat,
  SE = se,
  CI_low = ci[1],
  CI_high = ci[2]
)
print(results)



# Run again without covariates
sl_lib0 <- list("SL.mean")
X_empty <- pension[, c()]   #covariates


Q0 <- rep(NA, N)
g <- rep(NA, N)

coefs.Q0 <- rep(0, length(sl_lib0))
coefs.g <- rep(0, length(sl_lib0))

for(v in 1:nsplits) {
  train <- split_id != v
  test  <- split_id == v
  
  Q0.sl <- SuperLearner(Y[(A == 0) & train], X_empty[(A == 0) & train, , drop=FALSE], newX = X_empty[test, , drop=FALSE], family = gaussian(), SL.library = sl_lib0)
  Q0[test] <- Q0.sl$SL.predict
  coefs.Q0 <- coefs.Q0 + Q0.sl$coef
  
  g.sl <- SuperLearner(A[train], X_empty[train, , drop=FALSE], newX = X_empty[test, , drop=FALSE], family= binomial(), SL.library = sl_lib0)
  g[test] <- g.sl$SL.predict
  coefs.g <- coefs.g + g.sl$coef
}


P_A1 <- mean(A)
g <- pmin(g, 0.99)
g <- pmax(g, 0.01)

# Construct clever covariate H
H <- ifelse(A == 1, 1 / P_A1, -g / (P_A1 * (1 - g)))

# Fit fluctuation model (update Q0 to Q0_star)
epsilon_model <- glm(Y ~ -1 + H + offset(Q0), family = gaussian())
epsilon <- coef(epsilon_model)["H"]
Q0_tmle <- Q0 + epsilon * H

# Compute ATT
theta_hat <- mean( (Y[A == 1] - Q0_tmle[A == 1]) )

# Influence function and standard error
IF_part1 <- (A / P_A1) * (Y - Q0_tmle - theta_hat)
IF_part2 <- ((1 - A) * g / (P_A1 * (1 - g))) * (Y - Q0_tmle)
IF <- IF_part1 - IF_part2
se <- sqrt(mean(IF^2)) / sqrt(N)

# Confidence interval
ci <- theta_hat + c(-1.96, 1.96) * se

# Results
results <- data.frame(
  ATT = theta_hat,
  SE = se,
  CI_low = ci[1],
  CI_high = ci[2]
)
print(results)

"""
results:
       ATT       SE   CI_low  CI_high
1 7266.492 1416.717 4489.727 10043.26
"""



### QUESTION 3 (Florin)

## Prep
Y <- pension$net_tfa
A <- pension$e401
L <- subset(pension, select = c(age, inc, fsize, educ, marr, twoearn, db, pira, hown))

N <- dim(pension)[1]

SL.library <- c("SL.glm", "SL.step", "SL.ranger")

## Influence curve (IC), estimator, SE and 95% CI calculations using 5-fold cross-fitting
ICn = vector("numeric",length=N)

for (k in 1:5){
  evalset = (1+(N/5)*(k-1)):(N*k/5)
  
  p.model = SuperLearner(Y = A[-evalset], X = L[-evalset,], SL.library = SL.library, family = binomial)
  p = predict(p.model,newdata = L[evalset,])$pred
  q.model = SuperLearner(Y = Y[-evalset], X = L[-evalset,], SL.library = SL.library)
  q = predict(q.model, newdata = L[evalset,])$pred
  
  ICn[evalset] = (A[evalset] - p) * (Y[evalset] - q)
}

theta = mean(ICn)
IC = ICn-theta
SE = sd(IC)/sqrt(N)

theta_lb <- theta - 1.96 * SE
theta_ub <- theta + 1.96 * SE

## Testing null hypothesis that theta is 0
t.statistic <- t.test(x = ICn, alternative = c("two.sided"), mu = 0)$statistic
df <- t.test(x = ICn, alternative = c("two.sided"), mu = 0)$parameter
p.value <- t.test(x = ICn, alternative = c("two.sided"), mu = 0)$p.value

## Combining results into a table
results <- tibble::tibble(
  Method = "Eff. Influence Function",
  Estimate = theta,
  StdError = SE,
  CI_Lower = theta_lb,
  CI_Upper = theta_ub,
  T_value = t.statistic,
  P_value = p.value,
  Df = df
)

print(results)

## Quality checks
hist(ICn)




### QUESTION 4 (Thomas)

"""
Use the R-learner to estimate the effect of eligibility for enrolling in a 401(k) plan on net financial assets 
in function of age, income and years of education, but controlling for the variables listed in question 2. 
Report a histogram of the effect estimates obtained for all individuals in the study. 
In addition, visualize how the effects vary by age, income and years of education.
"""
#starting from lab 7 exercise 5


#initialize
#SL.library <- c("SL.glm","SL.glm.interaction","SL.ranger")


SL.library <- c("SL.glm", "SL.ranger", "SL.glmnet")


A <- pension$e401      #eligibility (treatment)
Y <- pension$net_tfa   #net financial assets (outcome)
L <- pension[, c("age", "inc", "fsize", "educ", "marr", "twoearn", "db", "pira", "hown" )]   #covariates
Z <- pension[, c("age","inc","educ")]

#Train a prediction model

# Create 5 folds
set.seed(123)
N <- length(A) #N = 9915
K <- 5
#each unit is randomly assigned to fold (from 1 to 5) 
fold_id <- sample(rep(1:K, length.out = N))  # = 9915

m_hat <- numeric(N)
e_hat <- numeric(N)


for(k in 1:K){
  #create training and validation set
  train_idx <- which(fold_id != k)
  valid_idx <- which(fold_id == k)
  
  fit_m = SuperLearner(Y = Y[train_idx], X = L[train_idx, , drop = FALSE], SL.library = SL.library)
  fit_e =  SuperLearner(Y = A[train_idx], X = L[train_idx, , drop = FALSE], SL.library = SL.library, family = binomial())
  
  #predict on validation fold
  newL <- L[valid_idx, , drop = FALSE]
  m_hat[valid_idx] = predict(fit_m, newdata = newL)$pred
  e_hat[valid_idx] = predict(fit_e, newdata = newL)$pred
}


# Pseudo-outcome and weights
pseudoY <- (Y - m_hat) / (A - e_hat)
W      <- (A - e_hat)^2


# Predict in function of age, income and years of education
fit_tau <- SuperLearner(
  Y          = pseudoY,
  X          = Z,
  SL.library = SL.library,
  obsWeights = W
)
tau_hat <- predict(fit_tau, newdata = Z)$pred


## estimate + CI
N         <- length(tau_hat)
mean_tau  <- mean(tau_hat)
se_tau    <- sd(tau_hat) / sqrt(N)

# 95% CI via normal approximation
ci_tau    <- mean_tau + c(-1.96, +1.96) * se_tau

# Or even more simply:
t_res     <- t.test(tau_hat)
mean_tau  <- t_res$estimate
ci_tau    <- t_res$conf.int

# Print
cat(sprintf("Average CATE = %.2f (95%% CI: [%.2f, %.2f])\n",
            mean_tau, ci_tau[1], ci_tau[2]))



#check for outliers
summary(tau_hat)
quantile(tau_hat, probs = c(0.005, 0.025, 0.975, 0.995))

summary(e_hat)
hist(e_hat, breaks=20, main="Propensity score distribution")

lims <- quantile(tau_hat, c(0.01, 0.99))


# Histogram of individual CATEs:
hist(tau_hat,
     xlab = "Estimated CATE (R-learner)",
     main = "401(k) effects")


# Clipped version of previous plot
qs     <- quantile(tau_hat, c(0.005, 0.995))
tau_clipped <- pmin(pmax(tau_hat, qs[1]), qs[2])
hist(
  tau_clipped, 
  breaks = 50, 
  xlab   = "Winsorized CATE",
  main   = "401(k) effects (0.5–99.5% winsorized)"
)



SL.library.visual <- c("SL.glm", "SL.ranger", "SL.gam")

# visualization of how the effects vary by age, income and years of education.
fit_age <- SuperLearner(Y = pseudoY, X = Z[,"age",  drop=FALSE], SL.library = SL.library.visual, obsWeights = W)

grid_age <- data.frame(age = seq(min(L$age), max(L$age), length=100))
tau_age <- predict(fit_age, newdata = grid_age, onlySL = TRUE)$pred
plot(grid_age$age, tau_age, type="l",
     xlab="Age", ylab="E[τ | age]",
     main="CATE Age")

plot(
  grid_age$age, tau_age, 
  type  = "l",
  lwd   = 2,
  xlab  = "Age",
  ylab  = expression(hat(tau)(age)),
  main  = "CATE Age"
)


#for income
fit_inc <- SuperLearner(Y = pseudoY, X = Z[,"inc",  drop=FALSE], SL.library = SL.library.visual, obsWeights = W)

grid_inc <- data.frame(inc = seq(min(Z$inc), max(Z$inc), length=100))
tau_inc <- predict(fit_inc, newdata = grid_inc, onlySL = TRUE)$pred
plot(grid_inc$inc, tau_inc, type="l",
     xlab="Income", ylab="E[τ | income]",
     main="Marginal CATE vs Income")

plot(
  grid_inc$inc, tau_inc, 
  type  = "l",
  lwd   = 2,
  xlab  = "Income",
  ylab  = expression(hat(tau)(inc)),
  main  = "CATE Income"
)


#for years of education
fit_educ <- SuperLearner(Y = pseudoY, X = Z[,"educ",  drop=FALSE], SL.library = SL.library.visual, obsWeights = W)

grid_educ <- data.frame(educ = seq(min(Z$educ), max(Z$educ), length=100))
tau_educ <- predict(fit_educ, newdata = grid_educ, onlySL = TRUE)$pred
plot(grid_educ$educ, tau_educ, type="l",
     xlab="Education", ylab="E[τ | educ]",
     main="CATE Education")

plot(
  grid_educ$educ, tau_educ, 
  type  = "l",
  lwd   = 2,
  xlab  = "Years of education",
  ylab  = expression(hat(tau)(educ)),
  main  = "CATE Education"
)

min(pension$inc)
table(L$educ)[as.character(0:18)]


### QUESTION 5 (Florin + Thomas)

## Prep
Z <- subset(L, select = c(age, inc, educ))

## DR-learner predictions using 5-fold cross-fitting
Cn = vector("numeric",length=N)
Ytilde1 = vector("numeric",length=N)

for (k in 1:5){
  evalset = (1+(N/5)*(k-1)):(N*k/5)
  Y.nonevalset = Y[-evalset]
  L.nonevalset = L[-evalset,]
  A.nonevalset = A[-evalset]
  
  p.model = SuperLearner(Y = A.nonevalset, X = L.nonevalset, SL.library = SL.library, family = binomial)
  p = predict(p.model, newdata = L[evalset,])$pred
  q1.model = SuperLearner(Y = Y.nonevalset[A.nonevalset == 1], X = L.nonevalset[A.nonevalset == 1,], SL.library = SL.library)
  q1 = predict(q1.model, newdata = L[evalset,])$pred
  
  Cn[evalset] = q1 + A[evalset]/p * (Y[evalset] - q1)
  pseudo.outcome.model = SuperLearner(Y = Cn[evalset], X = Z[evalset,], SL.library = SL.library)
  Ytilde1[evalset] = predict(pseudo.outcome.model, newdata = Z[evalset,])$pred
}

hist(Ytilde1, xlab = "DR-learner predictions", main = "")

## Estimate of counterfactual prediction error with 95% confidence interval (5-fold cross-fitting)
ICn = vector("numeric",length=N)

for (k in 1:5){
  evalset = (1+(N/5)*(k-1)):(N*k/5)
  Y.nonevalset = Y[-evalset]
  L.nonevalset = L[-evalset,]
  A.nonevalset = A[-evalset]
  
  Ytilde1.nonevalset = Ytilde1[-evalset]
  
  p.model = SuperLearner(Y = A.nonevalset, X = L.nonevalset, SL.library = SL.library, family = binomial)
  p = predict(p.model, newdata = L[evalset,])$pred
  q1.model = SuperLearner(Y = Y.nonevalset[A.nonevalset == 1], X = L.nonevalset[A.nonevalset == 1,], SL.library = SL.library)
  q1 = predict(q1.model, newdata = L[evalset,])$pred
  
  qtilde1.model = SuperLearner(Y = Ytilde1.nonevalset[A.nonevalset == 1], X = L.nonevalset[A.nonevalset == 1,], SL.library = SL.library)
  qtilde1 = predict(q1.model, newdata = L[evalset,])$pred
  
  ICn[evalset] = ((A[evalset]/p) * (Y[evalset] - q1) + q1)^2 + ((A[evalset]/p) * (Ytilde1[evalset] - qtilde1) + qtilde1)^2 # slightly adjusted influence function for ATE to (hopefully) estimate the counterfactual prediction error in a debiased manner 
}

est = mean(ICn)
IC = ICn-est
SE = sd(IC)/sqrt(N)

est_lb <- est - 1.96 * SE
est_ub <- est + 1.96 * SE

## Combining results into a table
results <- tibble::tibble(
  Method = "Eff. Influence Function",
  Estimate = est,
  StdError = SE,
  CI_Lower = est_lb,
  CI_Upper = est_ub
)

results <-
  results %>%
  mutate(across(where(is.numeric), ~ num(., digits = 3)))

print(results)

