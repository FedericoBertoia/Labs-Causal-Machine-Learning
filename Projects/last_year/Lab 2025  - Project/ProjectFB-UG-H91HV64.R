###################
# Setup
###################

library(hdm)
data(pension)

set.seed(123)

#########################################################
# Question 2 - AIPW and TMLE for ATE and ATT
#########################################################

library(SuperLearner)
library(caret)

# - - - Variables - - - 

Y <- pension$net_tfa
A <- pension$e401
L <- pension[,c("age", "inc", "fsize", "educ", "marr", "twoearn", "db", "ira", "hown")]
X <- cbind(A, L)

N <- dim(pension)[1]
folds <- createFolds(A, k=5)
SL.library <- c("SL.glm", "SL.earth", "SL.ranger")

# - - - Nuisance Parameters Estimation - - -

Y1 <- rep(NA, N)
Y0 <- rep(NA, N)
pi <- rep(NA, N)

for (k in 1:5){
  test <- folds[[k]]
  train <- setdiff(1:N, test)
  
  model.Y <- SuperLearner(Y[train], X[train,], SL.library = SL.library, family = gaussian())
  
  model.A <- SuperLearner(A[train], L[train,], SL.library = SL.library, family = binomial())
  
  Y1[test] <- predict(model.Y, newdata = data.frame(A = 1, L[test,]), onlySL = TRUE)$pred
  Y0[test] <- predict(model.Y, newdata = data.frame(A = 0, L[test,]), onlySL = TRUE)$pred
  
  pi[test] <- predict(model.A, newdata = data.frame(L[test,]), onlySL = TRUE)$pred
  
}

Q <- A*Y1 + (1-A)*Y0

# - - - Check stability of inverse propensity scores - - -
pi_A1 <- pi[A==1]
pi_A0 <- pi[A==0]

pi_A1_inverse <- 1/pi_A1
pi_A0_inverse <- 1/pi_A0

# - - - Distribution propensity scores - - -

df_density <- data.frame(
  pi = c(pi_A1, pi_A0),
  A = factor(rep(c(1, 0), c(length(pi_A1), length(pi_A0))))
)

ggplot(df_density, aes(x = pi, color = A, fill = A)) +
  geom_density(alpha = 0.3) +
  labs(
    title = "Kernel Density of Propensity Scores",
    x = "Propensity Score", y = "Density",
    color = "Treatment A", fill = "Treatment A"
  ) +
  scale_color_manual(values = c("red", "blue")) +
  scale_fill_manual(values = c("red", "blue")) +
  theme_minimal(base_size = 14)



# - - - Boxplot of propensity scores - - -

df_box_pi <- data.frame(
  pi = c(pi_A1, pi_A0),
  A = factor(rep(c(1, 0), c(length(pi_A1), length(pi_A0))))
)

ggplot(df_box_pi, aes(x = A, y = pi, fill = A)) +
  geom_boxplot() +
  labs(
    title = "Boxplot of Propensity Scores",
    x = "Treatment A", y = "Propensity Score"
  ) +
  scale_fill_manual(values = c("red", "blue")) +
  theme_minimal(base_size = 14)



# - - - Boxplot of (inverse) propensity scores - - -

df_box_inv <- data.frame(
  inv_pi = c(pi_A1_inverse, pi_A0_inverse),
  A = factor(rep(c(1, 0), c(length(pi_A1_inverse), length(pi_A0_inverse))))
)

ggplot(df_box_inv, aes(x = A, y = inv_pi, fill = A)) +
  geom_boxplot() +
  labs(
    title = "Boxplot of Inverse Propensity Scores",
    x = "Treatment A", y = "1 / Propensity Score"
  ) +
  scale_fill_manual(values = c("red", "blue")) +
  theme_minimal(base_size = 14)


# - - - Round propensity scores - - -

#pi <- ifelse(pi < 0.05, 0.05, ifelse(pi > 0.95, 0.95, pi))



# - - - AIPW - - -

# ATE

ate_aipw <- mean(((A/pi) - ((1-A)/(1-pi))) * (Y - Q) + Y1 - Y0)
ate_aipw_eic <- ((A/pi) - ((1-A)/(1-pi))) * (Y - Q) + Y1 - Y0 - ate_aipw

ate_aipw_se <- (1/sqrt(N)) * sd(ate_aipw_eic)

ate_aipw_lb <- ate_aipw - 1.96 * ate_aipw_se
ate_aipw_ub <- ate_aipw + 1.96 * ate_aipw_se

boxplot(ate_aipw_eic,
        main = "Zoomed Boxplot of AIPW EIC ATE",
        ylab = "EIC values",
        col = "lightgreen",
        ylim = quantile(ate_aipw_eic, c(0.05, 0.95)))



# ATT

p_A <- mean(A)

att_aipw <- mean(((A/p_A) * (Y - Y0)) - (((1-A)/(1-pi)) * (pi/p_A) * (Y - Y0)))

att_aipw_eic <- ((A/p_A) * (Y - Y0 - att_aipw)) - (((1-A)/(1-pi)) * (pi/p_A) * (Y - Y0))

att_aipw_se <- (1/sqrt(N)) * sd(att_aipw_eic)

att_aipw_lb <- att_aipw - 1.96 * att_aipw_se
att_aipw_ub <- att_aipw + 1.96 * att_aipw_se


boxplot(att_aipw_eic,
        main = "Zoomed Boxplot of AIPW EIC ATT",
        ylab = "EIC values",
        col = "lightgreen",
        ylim = quantile(att_aipw_eic, c(0.05, 0.95)))

# - - - TMLE - - -

library(tmle)

tmle <- tmle(Y = Y, A = A, W = L, Q = cbind(Y0, Y1), g1W = pi, family = "gaussian")

summary(tmle)

# - - - TMLE BY HAND (CORRECTED) - - -

# Fluctuation model: linear regression for control units only (A=0)
fluctmod <- glm(Y[A == 0] ~ -1 + offset(Y0[A == 0]) + I(pi[A == 0] / (1 - pi[A == 0])), family = gaussian)
epsilon <- coef(fluctmod)

# Update Q0 for all units using fluctuation coefficient times clever covariate
Q0_tmle <- Y0 + epsilon * (pi / (1 - pi))

# Compute ATT estimator (weighted mean difference for treated)
theta_hat <- mean((A / p_A) * (Y - Q0_tmle))

# Influence curve (IC)
IC <- (A / p_A) * (Y - Q0_tmle - theta_hat) - ((1 - A) * pi / (1 - pi)) * (Y - Q0_tmle) / p_A

# Standard error
se <- sd(IC) / sqrt(N)

# Confidence interval
ci <- theta_hat + c(-1.96, 1.96) * se

# Results
list(ATT = theta_hat, SE = se, CI = ci)



# - - - Unadjusted Analysis - - -

unadjusted = net_tfa ~ e401 

model = lm(unadjusted, data = pension)
summary(model)
confint(model)


# --- Combine results into a table ---                                             
results_ate <- tibble::tibble(
  Method = c("AIPW", "TMLE", "Unadjusted"),
  Estimate = c(ate_aipw, tmle$estimates$ATE$psi, model$coefficients[2]),
  StdError = c(ate_aipw_se, sqrt(tmle$estimates$ATE$var.psi), 1305.7),
  CI_Lower = c(ate_aipw_lb, tmle$estimates$ATE$CI[1], 16999.9),
  CI_Upper = c(ate_aipw_ub, tmle$estimates$ATE$CI[2], 22118.8)
)

print(results_ate)


# --- Combine results into a table ---                                             
results_att <- tibble::tibble(
  Method = c("AIPW", "TMLE"),
  Estimate = c(att_aipw, tmle$estimates$ATT$psi),
  StdError = c(att_aipw_se, sqrt(tmle$estimates$ATT$var.psi)),
  CI_Lower = c(att_aipw_lb, tmle$estimates$ATT$CI[1]),
  CI_Upper = c(att_aipw_ub, tmle$estimates$ATT$CI[2])
)

print(results_att)

#########################################################
# Question 3 - Test expected covariance = 0 using EIF
#########################################################
Y_hat <- rep(NA,N)

for (k in 1:5){
  test <- folds[[k]]
  train <- setdiff(1:N, test)
  
  model.YL <- SuperLearner(Y[train], L[train,], SL.library = SL.library, family = gaussian())
  
  Y_hat[test] <- predict(model.YL, newdata = data.frame(L[test,]), onlySL = TRUE)$pred

}

# PROVA PER VEDERE COS'HANNO FATTO ALL'ESAME
mu <- (A*Y1) + (1-A)*(Y0)

theta2 <- mean((A - pi)*(Y - mu)); theta2
# PROVA PER VEDERE COS'HANNO FATTO ALL'ESAME

theta <- mean((A - pi)*(Y - Y_hat));theta

EIC_theta <- (A-pi)*(Y-Y_hat) - theta

theta_se <- (1/sqrt(N)) * sd(EIC_theta)

# 95% CI
theta_ci_lb <- theta - 1.96 * theta_se
theta_ci_ub <- theta + 1.96 * theta_se

boxplot(EIC_theta,
        main = "Zoomed Boxplot of EIC_theta",
        ylab = "EIC values",
        col = "lightgreen",
        ylim = quantile(EIC_theta, c(0.05, 0.95)))


# --- Combine results into a table ---
results3 <- tibble::tibble(
  Method = c("Eff. Influence Function"),
  Estimate = c(theta),
  StdError = c(theta_se),
  CI_Lower = c(theta_ci_lb),
  CI_Upper = c(theta_ci_ub)
)

print(results3)

#########################################################
# Question 4 - R Learner for CATE
#########################################################

L2 <- pension[,c("age", "inc", "educ")]

# --- Approach 1 (No cross fitting for last step) ---

Y_res <- Y - Y_hat
W_res <- A - pi

pseudo_R <- Y_res/W_res
weights <- (W_res)^2

model.R_learner <- SuperLearner(Y = pseudo_R, X = L2, SL.library = SL.library, family = gaussian(), obsWeights = weights)

cate.R_learner1 <- predict(model.R_learner, newdata = data.frame(L2), onlySL = TRUE)$pred

summary(cate.R_learner1)

cate_core1 <- cate.R_learner1[cate.R_learner1 < 30000 & cate.R_learner1 > -10000]

hist(cate_core1,
     breaks = 100,
     col = "lightgreen",
     border = "white",
     main = "Histogram of CATE (R-learner 1, trimmed)",
     xlab = "CATE Estimates")

# AGE

model.R_learner_age <- SuperLearner(Y = pseudo_R, X = as.data.frame(pension$age), SL.library = SL.library, family = gaussian(), obsWeights = weights)

cate.R_learner1_age <- predict(model.R_learner_age, newdata = as.data.frame(pension$age), onlySL = TRUE)$pred

summary(cate.R_learner1_age)

library(ggplot2)

df <- data.frame(
  age = pension$age,
  cate = cate.R_learner1_age
)

ggplot(df, aes(x = age, y = cate)) +
  geom_point(alpha = 0.5, color = "steelblue") +
  labs(
    title = "CATE vs Age",
    x = "Age",
    y = "CATE (R-learner)"
  ) +
  theme_minimal()




# INCOME

model.R_learner_income <- SuperLearner(Y = pseudo_R, X = as.data.frame(pension$inc), SL.library = SL.library, family = gaussian(), obsWeights = weights)

cate.R_learner1_income <- predict(model.R_learner_income, newdata = as.data.frame(pension$inc), onlySL = TRUE)$pred

summary(cate.R_learner1_income)

df <- data.frame(
  age = pension$inc,
  cate = cate.R_learner1_income
)

ggplot(df, aes(x = age, y = cate)) +
  geom_point(alpha = 0.5, color = "steelblue") +
  labs(
    title = "CATE vs Income",
    x = "Income",
    y = "CATE (R-learner)"
  ) +
  theme_minimal()



# EDUCATION

model.R_learner_educ <- SuperLearner(Y = pseudo_R, X = as.data.frame(pension$educ), SL.library = SL.library, family = gaussian(), obsWeights = weights)

cate.R_learner1_educ <- predict(model.R_learner_educ, newdata = as.data.frame(pension$educ), onlySL = TRUE)$pred

summary(cate.R_learner1_educ)

df <- data.frame(
  age = pension$educ,
  cate = cate.R_learner1_educ
)

ggplot(df, aes(x = age, y = cate)) +
  geom_point(alpha = 0.5, color = "steelblue") +
  labs(
    title = "CATE vs Education",
    x = "Education",
    y = "CATE (R-learner)"
  ) +
  theme_minimal()




# --- Approach 2 (Cross fitting for last step as in question 5) ---

folds3 <- createFolds(A, k=3)
cate.R_learner2 <- rep(NA, N)
cate.R_learner2_age <- rep(NA, N)
cate.R_learner2_inc <- rep(NA, N)
cate.R_learner2_educ <- rep(NA, N)

for (k in 1:3){
  test <- folds3[[k]]
  others <- setdiff(1:3, k)
  trainA <- folds3[[others[1]]]
  trainB <- folds3[[others[2]]]
  
  # Nuisance Fit
  model.YL.trainA <- SuperLearner(Y=Y[trainA], X <- L[trainA,], SL.library = SL.library, family = gaussian())
  
  model.YL.trainB <- SuperLearner(Y=Y[trainB], X <- L[trainA,], SL.library = SL.library, family = gaussian())
  
  model.pi.trainA <- SuperLearner(Y=A[trainA], X <- L[trainA,], SL.library = SL.library, family = binomial())
  model.pi.trainB <- SuperLearner(Y=A[trainB], X <- L[trainB,], SL.library = SL.library, family = binomial())
  
  # Nuisance Pred
  Y_hat.trainA <- predict(model.YL.trainB, newdata = L[trainA,], onlySL = TRUE)$pred
  
  Y_hat.trainB <- predict(model.YL.trainA, newdata = L[trainB,], onlySL = TRUE)$pred
  
  pi.trainA <- predict(model.pi.trainB, newdata = L[trainA,], onlySL = TRUE)$pred
  pi.trainB <- predict(model.pi.trainA, newdata = L[trainB,], onlySL = TRUE)$pred
  
  # Pseudo Outcome
  Y_res.trainA <- Y[trainA] - Y_hat.trainA
  Y_res.trainB <- Y[trainB] - Y_hat.trainB
  Y_res.trainAB <- c(Y_res.trainA, Y_res.trainB)
  
  W_res.trainA <- A[trainA] - pi.trainA
  W_res.trainB <- A[trainB] - pi.trainB
  W_res.trainAB <- c(W_res.trainA, W_res.trainB)
  
  pseudo_R.trainAB <- (Y_res.trainAB / W_res.trainAB)
  weights.trainAB <- (W_res.trainAB)^2
  
  #ALL
  model.pseudo <- SuperLearner(Y = pseudo_R.trainAB, X = as.data.frame(L2[c(trainA, trainB), ]), SL.library = SL.library, family = gaussian(), obsWeights = weights.trainAB)
  
  cate.R_learner2[test] <- predict(model.pseudo, newdata = as.data.frame(L2[test,]), onlySL = TRUE)$pred
  
  # AGE
  model.pseudo_age <- SuperLearner(Y = pseudo_R.trainAB, X = as.data.frame(L2[c(trainA, trainB), "age"]), SL.library = SL.library, family = gaussian(), obsWeights = weights.trainAB)
  
  cate.R_learner2_inc[test] <- predict(model.pseudo_age, newdata = as.data.frame(L2[test, "age"]), onlySL = TRUE)$pred
  
  # INCOME
  model.pseudo_inc <- SuperLearner(Y = pseudo_R.trainAB, X = as.data.frame(L2[c(trainA, trainB), "inc"]), SL.library = SL.library, family = gaussian(), obsWeights = weights.trainAB)
  
  cate.R_learner2_inc[test] <- predict(model.pseudo_inc, newdata = as.data.frame(L2[test, "inc"]), onlySL = TRUE)$pred
  
  # EDUC
  model.pseudo_educ <- SuperLearner(Y = pseudo_R.trainAB, X = as.data.frame(L2[c(trainA, trainB), "educ"]), SL.library = SL.library, family = gaussian(), obsWeights = weights.trainAB)
  
  cate.R_learner2_educ[test] <- predict(model.pseudo_educ, newdata = as.data.frame(L2[test, "educ"]), onlySL = TRUE)$pred
  
}

# Problem because X is not a matrix, fix later...



#########################################################
# Question 5 - DR Learner for counterfactual error
#########################################################

folds3 <- createFolds(A, k=3)              
Y1_tilde <- rep(NA, N)

for (k in 1:3){
  test <- folds3[[k]]
  others <- setdiff(1:3, k)
  trainA <- folds3[[others[1]]]
  trainB <- folds3[[others[2]]]
  
  # Nuisance Fit
  model.Y1.trainA <- SuperLearner(Y=Y[trainA][A[trainA] ==1], X <- L[trainA,][A[trainA] == 1,], SL.library = SL.library, family = gaussian())
  model.Y0.trainA <- SuperLearner(Y=Y[trainA][A[trainA] ==0], X <- L[trainA,][A[trainA] == 0,], SL.library = SL.library, family = gaussian())
  
  model.Y1.trainB <- SuperLearner(Y=Y[trainB][A[trainB] ==1], X <- L[trainA,][A[trainB] == 1,], SL.library = SL.library, family = gaussian())
  model.Y0.trainB <- SuperLearner(Y=Y[trainB][A[trainB] ==0], X <- L[trainA,][A[trainB] == 0,], SL.library = SL.library, family = gaussian())
  
  model.pi.trainA <- SuperLearner(Y=A[trainA], X <- L[trainA,], SL.library = SL.library, family = binomial())
  model.pi.trainB <- SuperLearner(Y=A[trainB], X <- L[trainB,], SL.library = SL.library, family = binomial())
  
  # Nuisance Pred
  Y1.trainA <- predict(model.Y1.trainB, newdata = L[trainA,], onlySL = TRUE)$pred
  Y0.trainA <- predict(model.Y0.trainB, newdata = L[trainA,], onlySL = TRUE)$pred
  
  Y1.trainB <- predict(model.Y1.trainA, newdata = L[trainB,], onlySL = TRUE)$pred
  Y1.trainB <- predict(model.Y0.trainA, newdata = L[trainB,], onlySL = TRUE)$pred
  
  pi.trainA <- predict(model.pi.trainB, newdata = L[trainA,], onlySL = TRUE)$pred
  pi.trainB <- predict(model.pi.trainA, newdata = L[trainB,], onlySL = TRUE)$pred
  
  # Pseudo Outcome
  
  Ci.trainA <- Y1.trainA + (A[trainA]/pi.trainA)*(Y[trainA]-Y1.trainA)
  Ci.trainB <- Y1.trainB + (A[trainB]/pi.trainB)*(Y[trainB]-Y1.trainB)
  Ci.trainAB <- c(Ci.trainA, Ci.trainB)
  
  model.pseudo <- SuperLearner(Y = Ci.trainAB, X = L[c(trainA, trainB), ], SL.library = SL.library, family = gaussian())
  
  Y1_tilde[test] <- predict(model.pseudo, newdata = L[test,], onlySL = TRUE)$pred
  
}

summary(Y1_tilde)

Y1_tilde_core <- Y1_tilde[Y1_tilde < 50000 & Y1_tilde > -25000]

hist(Y1_tilde_core,
     breaks = 100,
     col = "lightgreen",
     border = "white",
     main = "Histogram of Y1_Tilde (DR-learner)",
     xlab = "Y1_Tilde Estimates")

#################
# Plug In
#################

# Estimate Prediction Error and 95% Confidence Interval
mse_est <- mean((Y[A==1] - Y1_tilde[A==1])^2)

# Standard error and 95% CI
se <- sd((Y[A==1] - Y1_tilde[A==1])^2) / sqrt(sum(A))
ci_low <- mse_est - 1.96 * se
ci_high <- mse_est + 1.96 * se

cat("Estimated prediction error:", mse_est, "\n")

cat("95% CI: [", ci_low, ",", ci_high, "]\n")

#################
# Based on the EIF
#################

Z = (Y - Y1_tilde)^2

Z_hat <- rep(NA, N)


for (k in 1:5){
  test <- folds[[k]]
  train <- setdiff(1:N, test)
  
  sl.Z <- SuperLearner(Y = Z[train][A[train]==1], X = L[train,][A[train]==1,], SL.library = SL.library, family = gaussian())
  
  z_pred <- predict(sl.Z, newdata = L[test,], onlySL = TRUE)$pred
  
  Z_hat[test] <- z_pred
} 


mse_eif <- mean(((A*(Z - Z_hat))/pi) + Z_hat)

eic <- ((A*(Z - Z_hat)/pi)) + Z_hat - mse_eif

se_mse_eif <- (1/sqrt(N)) * sd(eic)


mse_eif_lb <- mse_eif - 1.96*se_mse_eif
mse_eif_ub <- mse_eif + 1.96*se_mse_eif

# --- Combine results into a table ---
results <- tibble::tibble(
  Method = c("Plug-In", "Eff. Influence Function"),
  Estimate = c(mse_est, mse_eif),
  StdError = c(se, se_mse_eif),
  CI_Lower = c(ci_low, mse_eif_lb),
  CI_Upper = c(ci_high, mse_eif_ub)
)

print(results)


# PROVA PER VEDERE COS'HANNO FATTO ALL'ESAME
ciao <- mean((A * Z_hat) / pi)
# PROVA PER VEDERE COS'HANNO FATTO ALL'ESAME
