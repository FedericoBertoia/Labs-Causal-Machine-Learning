###############################################
### 	LAB 5: The Average Treatment Effect    ###
###############################################

library(twang)
library(tibble)
library(dplyr)

################
# Load data 
################
set.seed(123)

data(lindner)
?lindner
lindner$survival_binary<-ifelse(lindner$sixMonthSurvive==FALSE,1,0)


#######################################################
###   Question 1 - Risk Difference and Odds ratio	  ###
#######################################################

table(lindner$survival_binary,lindner$abcix)

unadjusted = survival_binary ~ abcix
adjusted = survival_binary ~ abcix + stent + height + female + diabetic + acutemi + ejecfrac + ves1proc

# Risk Difference
model <- lm(unadjusted, data = lindner)
summary(model)
confint(model)

# Odds Ratio
logit_model <- glm(unadjusted, data = lindner, family = 'binomial')
summary(logit_model)
confint(logit_model)

odds_ratio = exp(logit_model$coefficients[2]); odds_ratio


odds_ratio_LB = exp(confint(logit_model)[2,1]); odds_ratio_LB
odds_ratio_UB = exp(confint(logit_model)[2,2]); odds_ratio_UB



#######################################
###   Question 2 - G-Computation	  ###
#######################################

Y <- lindner$survival_binary
A <- lindner$abcix
X <- lindner[, c("stent", "height", "female", "diabetic", "acutemi", "ejecfrac", "ves1proc")]

dat <- as.data.frame(cbind(Y, A, X))
formula_gcomp <- as.formula(paste("Y ~ A +", paste(colnames(X), collapse = " + ")))

# G-Computation By Hand using parametric model (a)
fit <- glm(formula_gcomp, family = 'binomial', data = dat)

dat0 <- dat; dat0$A <- 0
dat1 <- dat; dat1$A <- 1

Q0 <- predict(fit, newdata = dat0, type = "response")
Q1 <- predict(fit, newdata = dat1, type = "response")

ate_gcomp_byhand <- mean(Q1 - Q0); ate_gcomp_byhand



# G-Computation using stdReg (b)
library(stdReg)

fit <- glm(formula_gcomp, family = 'binomial', data = dat)
fit.std <- stdGlm(fit = fit, data = dat, X = "A")
print(summary(fit.std))

ate_gcomp <- fit.std$est[2] - fit.std$est[1]; ate_gcomp
ate_gcomp_se <- sqrt(fit.std$vcov[1] + fit.std$vcov[4] - 2* fit.std$vcov[3]); ate_gcomp_se



# G-Computation Parametric with functions (b - alternative)

muhat_par <- estimate_mu_par(Y = Y,
                             A = A,
                             X = X,
                             family = 'binomial')

ate_gcomp_par <- est_gcomp(muhat = muhat_par); ate_gcomp_par
ate_gcomp_par_se <- bootstrap_se_gcomp(Y = Y, A = A, X = X, family = 'binomial', B = 200); ate_gcomp_par_se



# G-Computation Non Parametric with SuperLearner (d)
library(SuperLearner)
sl_lib <- list("SL.glm","SL.step","SL.glm.interaction","SL.randomForest")

AX <- as.data.frame(cbind(A, X))

Q.sl <- SuperLearner(Y = Y,
                     X = AX,
                     SL.library = sl_lib,
                     family = 'binomial')

Q1 <- predict(Q.sl, newdata = data.frame(A = 1, X))$pred
Q0 <- predict(Q.sl, newdata = data.frame(A = 0, X))$pred

ate_gcomp_sl <- mean(Q1 - Q0); ate_gcomp_sl
Q.sl$coef


# G-Computation Non Parametric with functions and cross-fitting (d - alternative)

muhat_est_nonpar <- estimate_mu_nonpar(Y = Y,
                                       A = A,
                                       X = X,
                                       sl_lib = sl_lib,
                                       nsplits = 5,
                                       family = 'binomial')

muhat_nonpar <- muhat_est_nonpar$preds
coefs_gcomp <- muhat_est_nonpar$coefs; coefs_gcomp

ate_gcomp_nonpar <- est_gcomp(muhat = muhat_nonpar); ate_gcomp_nonpar

# You can avoid running the bootstrap part, since it takes quite some time
# and we can not get valid standard error
ate_gcomp_nonpar_se <- bootstrap_se_gcomp(Y = Y, A = A, X = X, family = 'binomial',
                                         B = 50, nonpar = TRUE, sl_lib = sl_lib,
                                         nsplits = 3, stratified = FALSE); ate_gcomp_nonpar_se




#############################
###   Question 3 - IPW	  ###
#############################


# IPW using a parametric model and then a weighted linear regression (a)
library(sandwich)

AX <- as.data.frame(cbind(A, X))

formula_g <- as.formula(paste("A ~", paste(colnames(X), collapse = " + ")))
ghat <- glm(formula_g, family = 'binomial', AX)$fitted.values

ghat_bounded = bound_pi(ghat)

weights <- (A/ghat_bounded) + ((1-A)/(1-ghat_bounded))

fit_ipw <- lm(Y ~ A, weights = weights)

ate_ipw <- fit_ipw$coefficients[2]; ate_ipw
ate_ipw_se <- sqrt(vcovHC(fit_ipw, type = "HC")[2,2]); ate_ipw_se



# IPW parametric with functions (a - alternative)

pihat_par <- estimate_pi_par(A=A, X=X)
pihat_par_bounded <- bound_pi(pihat_par)

fit_ipw_par <- est_ipw(Y = Y,
                       A = A,
                       pi = pihat_par_bounded)

ate_ipw_par <- fit_ipw_par[1]; ate_ipw_par 
ate_ipw_par_sd <- fit_ipw_par[2]; ate_ipw_par_sd 




# IPW using a SuperLearner for propensity scores and then a weighted linear regression (b)
g.sl <- SuperLearner(Y = A,
                     X = X,
                     SL.library = sl_lib,
                     family = 'binomial')

ghat_sl <- predict(g.sl, newdata = X)$pred

ghat_sl_bounded <- bound_pi(ghat_sl)

weights_sl <- (A/ghat_sl_bounded) + ((1-A)/(1-ghat_sl_bounded))


fit_ipw_sl <- lm(Y ~ A, weights = weights_sl)

ate_ipw_nonpar <- fit_ipw_sl$coefficients[2]; ate_ipw_nonpar
ate_ipw_nonpar_se <- sqrt(vcovHC(fit_ipw_sl, type = "HC")[2,2]); ate_ipw_nonpar_se




# IPW nonparametric (SuperLearner) with functions (b - alternative)

pihat_nonpar_fit <- estimate_pi_nonpar(A = A,
                                        X = X,
                                        sl_lib = sl_lib,
                                        nsplits = 5)

pihat_nonpar <- pihat_nonpar_fit$preds

pihat_nonpar_bounded <- bound_pi(pihat_nonpar)

fit_ipw_nonpar <- est_ipw(Y = Y,
                          A = A,
                          pi = pihat_nonpar_bounded)

ate_ipw_nonpar <- fit_ipw_nonpar[1]; ate_ipw_nonpar 
ate_ipw_nonpar_se <- fit_ipw_nonpar[2]; ate_ipw_nonpar_se


# Check positivity assumption (c)
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




#############################
###   Question 4 - AIPW	  ###
#############################
N <- length(A)

# AIPW (a + b)
g1_aipw <- pihat_nonpar_bounded
g0_aipw <- 1-pihat_nonpar_bounded
Q1_aipw <- muhat_nonpar$mu1
Q0_aipw <- muhat_nonpar$mu0
Q_aipw <- muhat_nonpar$mu

ate_aipw_by_hand <- mean(((A/g1_aipw) - (1-A)/g0_aipw)*(Y - Q_aipw) + Q1_aipw - Q0_aipw)

influence_curve_aipw <- ((A/g1_aipw) - (1-A)/g0_aipw)*(Y - Q_aipw) + Q1_aipw - Q0_aipw - ate_aipw_by_hand

ate_aipw_by_hand_se <- (1/sqrt(N)) * sd(influence_curve_aipw)

ate_aipw_by_hand
ate_aipw_by_hand_se



# AIPW , using npcausal (c)

library(devtools) #install.packages("devtools")
library(npcausal) #install_github("ehkennedy/npcausal")

aipw_1cv <- ate(y = Y, a = A, x = X, nsplits = 1, sl.lib = sl_lib)

aipw_2cv <- ate(y = Y, a = A, x = X, nsplits = 2, sl.lib = sl_lib)

aipw_5cv <- ate(y = Y, a = A, x = X, nsplits = 5, sl.lib = sl_lib)



# AIPW, with functions (c - alternative)
ate_aipw_nonpar_est <- est_aipw(Y, A, pihat_nonpar_bounded, muhat_nonpar)

ate_aipw_nonpar <- ate_aipw_nonpar_est[1]; ate_aipw_nonpar
ate_aipw_nonpar_se <- ate_aipw_nonpar_est[2]; ate_aipw_nonpar_se






#############################
###   Question 5 - TMLE	  ###
#############################

# TMLE (a + b + c)
Q.sl <- SuperLearner(Y = Y,
                    X = AX,
                    SL.library = sl_lib,
                    family = 'binomial')

g.sl <- SuperLearner(Y = A,
                    X = X,
                    SL.library = sl_lib,
                    family = 'binomial')


Q1 <- predict(Q.sl, newdata = data.frame(A = 1, X))$pred
Q0 <- predict(Q.sl, newdata = data.frame(A = 0, X))$pred

Q <- A*Q1 + (1-A)*Q0

g <- predict(g.sl, newdata = X)$pred

H1 <- A/g
H0 <- (1-A)/(1-g)

H <- H1-H0

# Fluctuation Model 1

# Estimate parameters delta in the fluctuation models
delta1_tmle1 <- coef(glm(Y ~ -1 + offset(qlogis(Q1)) + H1, family = 'binomial'))
delta0_tmle1 <- coef(glm(Y ~ -1 + offset(qlogis(Q0)) + H0, family = 'binomial'))

# Targeting step
Q1_tmle1 <- plogis(qlogis(Q1) + delta1_tmle1/g)
Q0_tmle1 <- plogis(qlogis(Q0) + delta0_tmle1/(1-g))

# Compute E[Y^1] and E[Y^0]
EY1_tmle1 <- mean(Q1_tmle1); EY1_tmle1
EY0_tmle1 <- mean(Q0_tmle1); EY0_tmle1
ATE_tmle1 <- EY1_tmle1 - EY0_tmle1; ATE_tmle1

# Compute standard errors of E[Y^1] and E[Y^0] using the efficient influence curves
EIC_EY1_tmle1 <- (A / g) * (Y - Q1_tmle1) + Q1_tmle1 - EY1_tmle1
EIC_EY0_tmle1 <- ((1-A) / (1-g)) * (Y - Q0_tmle1) + Q0_tmle1 - EY0_tmle1

se_EY1_tmle1 <- (1/sqrt(N)) * sd(EIC_EY1_tmle1); se_EY1_tmle1
se_EY0_tmle1 <- (1/sqrt(N)) * sd(EIC_EY0_tmle1); se_EY0_tmle1
se_ATE_tmle1 <- (1/sqrt(N)) * sd(EIC_EY1_tmle1-EIC_EY0_tmle1); se_ATE_tmle1


# Compute the confidence intervals for E[Y^1] and E[Y^0]
EY1_tmle1_CI_low <- EY1_tmle1 - 1.96 * se_EY1_tmle1; print(EY1_tmle1_CI_low)
EY1_tmle1_CI_high <- EY1_tmle1 + 1.96 * se_EY1_tmle1; print(EY1_tmle1_CI_high)

EY0_tmle1_CI_low <- EY0_tmle1 - 1.96 * se_EY0_tmle1; print(EY0_tmle1_CI_low)
EY0_tmle1_CI_high <- EY0_tmle1 + 1.96 * se_EY0_tmle1; print(EY0_tmle1_CI_high)

ATE_tmle1_CI_low <- ATE_tmle1 - 1.96 * se_ATE_tmle1; print(ATE_tmle1_CI_low)
ATE_tmle1_CI_high <- ATE_tmle1 + 1.96 * se_ATE_tmle1; print(ATE_tmle1_CI_high)



# Fluctuation Model 2

# Estimate parameters delta in the fluctuation models
delta1_tmle2 <- coef(glm(Y ~ offset(qlogis(Q1)), family = 'binomial', weights = H1))
delta0_tmle2 <- coef(glm(Y ~ offset(qlogis(Q0)), family = 'binomial', weights = H0))

# Targeting step
Q1_tmle2 <- plogis(qlogis(Q1) + delta1_tmle2/g)
Q0_tmle2 <- plogis(qlogis(Q0) + delta0_tmle2/(1-g))

# Compute E[Y^1] and E[Y^0]
EY1_tmle2 <- mean(Q1_tmle2); EY1_tmle2
EY0_tmle2 <- mean(Q0_tmle2); EY0_tmle2
ATE_tmle2 <- EY1_tmle2 - EY0_tmle2; ATE_tmle2

# Compute standard errors of E[Y^1] and E[Y^0] using the efficient influence curves
EIC_EY1_tmle2 <- (A / g) * (Y - Q1_tmle2) + Q1_tmle2 - EY1_tmle2
EIC_EY0_tmle2 <- ((1-A) / (1-g)) * (Y - Q0_tmle2) + Q0_tmle2 - EY0_tmle2

se_EY1_tmle2 <- (1/sqrt(N)) * sd(EIC_EY1_tmle2); se_EY1_tmle2
se_EY0_tmle2 <- (1/sqrt(N)) * sd(EIC_EY0_tmle2); se_EY0_tmle2
se_ATE_tmle2 <- (1/sqrt(N)) * sd(EIC_EY1_tmle2-EIC_EY0_tmle2); se_ATE_tmle2

# Compute the confidence intervals for E[Y^1] and E[Y^0]
EY1_tmle2_CI_low <- EY1_tmle2 - 1.96 * se_EY1_tmle2; print(EY1_tmle2_CI_low)
EY1_tmle2_CI_high <- EY1_tmle2 + 1.96 * se_EY1_tmle2; print(EY1_tmle2_CI_high)

EY0_tmle2_CI_low <- EY0_tmle2 - 1.96 * se_EY0_tmle2; print(EY0_tmle2_CI_low)
EY0_tmle2_CI_high <- EY0_tmle2 + 1.96 * se_EY0_tmle2; print(EY0_tmle2_CI_high)

ATE_tmle2_CI_low <- ATE_tmle2 - 1.96 * se_ATE_tmle2; print(ATE_tmle2_CI_low)
ATE_tmle2_CI_high <- ATE_tmle2 + 1.96 * se_ATE_tmle2; print(ATE_tmle2_CI_high)









# TMLE using tmle package (d)
library(tmle)

tmle_package <- tmle(Y = Y,
                     A = A,
                     W = X,
                     Q.SL.library = sl_lib,
                     g.SL.library = sl_lib,
                     family = 'binomial',
                     V.Q = 5,
                     V.g = 5)

summary(tmle_package)


sl_lib2 <- list("glmnet", "SL.polymars", "SL.kernelKnn", "SL.glm")

tmle_package2 <- tmle(Y = Y,
                     A = A,
                     W = X,
                     Q.SL.library = sl_lib2,
                     g.SL.library = sl_lib2,
                     family = 'binomial',
                     V.Q = 10,
                     V.g = 10)

summary(tmle_package2)





# TMLE using functions (d - alternative)

ate_tmle_nonpar_est <- est_tmle(Y = Y,
                                A = A,
                                X = X,
                                muhat = muhat_nonpar,
                                pihat = pihat_nonpar)


ate_tmle_nonpar <- ate_tmle_nonpar_est[1]; ate_tmle_nonpar
ate_tmle_nonpar_se <- ate_tmle_nonpar_est[2]; ate_tmle_nonpar_se













# (Bonus) TMLE with sample splitting
N = length(A)
nsplits = 5
split_inds <- sample(rep(1:nsplits, ceiling(N/2))[1:N])

Q <- as_tibble(matrix(0, nrow = N, ncol = 2,
                          dimnames = list(c(), c("Q0", "Q1"))))
g <- rep(NA, N)

coefs.Q <- rep(0, length(sl_lib))
coefs.g <- rep(0, length(sl_lib))

for (vfold in 1:nsplits) {
  train <- split_inds != vfold
  test <- split_inds == vfold
  if (nsplits == 1) {
    train <- test
  }
  Q0.sl <- SuperLearner(Y[(A == 0) & train], X[(A == 0) & train, ],
                         newX = X[test, ], family = 'binomial', SL.library = sl_lib)
  
  Q[test, "Q0"] <- Q0.sl$SL.predict
  coefs.Q <- coefs.Q + Q0.sl$coef
  
  Q1.sl <- SuperLearner(Y[(A == 1) & train], X[(A == 1) & train, ],
                         newX = X[test, ], family = 'binomial', SL.library = sl_lib)
  
  Q[test, "Q1"] <- Q1.sl$SL.predict
  coefs.Q <- coefs.Q + Q1.sl$coef
  
  g.sl <- SuperLearner(A[train], X[train, ], newX = X[test, ],
                       family='binomial', SL.library = sl_lib)
  
  g[test] <- g.sl$SL.predict
  coefs.g <- coefs.g + g.sl$coef
}

coefs.Q <- coefs.Q/(2*nsplits); coefs.Q
coefs.g <- coefs.g/nsplits; coefs.g

H1 <- A/g
H0 <- (1-A)/(1-g)

# Estimate parameters delta in the fluctuation model and target the predictions

Q.tmle <- as_tibble(matrix(0, nrow = N, ncol = 2,
                      dimnames = list(c(), c("Q0.tmle", "Q1.tmle"))))


delta1.tmle <- coef(glm(Y ~ -1 + offset(qlogis(Q$Q1)) + H1, family = 'binomial'))
Q.tmle["Q1.tmle"]<- plogis(qlogis(Q$Q1) + delta1.tmle/g)

delta0.tmle <- coef(glm(Y ~ -1 + offset(qlogis(Q$Q0)) + H0, family = 'binomial'))
Q.tmle["Q0.tmle"]<- plogis(qlogis(Q$Q0) + delta0.tmle/(1-g))


mean(H1*(Y - Q.tmle$Q1.tmle))
mean(H0*(Y - Q.tmle$Q0.tmle))

# Compute E[Y^1] and E[Y^0]
EY1.tmle <- mean(Q.tmle$Q1.tmle); EY1.tmle
EY0.tmle <- mean(Q.tmle$Q0.tmle); EY0.tmle
ATE.tmle <- EY1.tmle - EY0.tmle; ATE.tmle


# Compute standard errors of E[Y^1] and E[Y^0] using the efficient influence curves
EIC.EY1.tmle <- (A / g) * (Y - Q.tmle$Q1.tmle) + Q.tmle$Q1.tmle - EY1.tmle
EIC.EY0.tmle <- ((1-A) / (1-g)) * (Y - Q.tmle$Q0.tmle) + Q.tmle$Q0.tmle - EY0.tmle

se.EY1.tmle <- (1/sqrt(N)) * sd(EIC.EY1.tmle); se.EY1.tmle
se.EY0.tmle <- (1/sqrt(N)) * sd(EIC.EY0.tmle); se.EY0.tmle
se.ATE.tmle <- (1/sqrt(N)) * sd(EIC.EY1.tmle-EIC.EY0.tmle); se.ATE.tmle


# Compute the confidence intervals for E[Y^1] and E[Y^0]
EY1.tmle.CI.low <- EY1.tmle - 1.96 * se.EY1.tmle; print(EY1.tmle.CI.low)
EY1.tmle.CI.high <- EY1.tmle + 1.96 * se.EY1.tmle; print(EY1.tmle.CI.high)

EY0.tmle.CI.low <- EY0.tmle - 1.96 * se.EY0.tmle; print(EY0.tmle.CI.low)
EY0.tmle.CI.high <- EY0.tmle + 1.96 * se.EY0.tmle; print(EY0.tmle.CI.high)

ATE.tmle.CI.low <- ATE.tmle - 1.96 * se.ATE.tmle; print(ATE.tmle.CI.low)
ATE.tmle.CI.high <- ATE.tmle + 1.96 * se.ATE.tmle; print(ATE.tmle.CI.high)






###################################################
###   Question Bonus - AIPW and TMLE for ATT	  ###
###################################################
N = length(A)
nsplits = 5
split_inds <- sample(rep(1:nsplits, ceiling(N/2))[1:N])


QA1 <- as_tibble(matrix(0, nrow = N, ncol = 2,
                           dimnames = list(c(), c("Q0", "Q1"))))
Q0 <- rep(NA, N)
g <- rep(NA, N)



coefs.Q0 <- rep(0, length(sl_lib))
coefs.g <- rep(0, length(sl_lib))

# Estimate Nuisance Parameters using Cross-Fitting
for (vfold in 1:nsplits) {
  train <- split_inds != vfold
  test <- split_inds == vfold
  if (nsplits == 1) {
    train <- test
  }
  Q0.sl <- SuperLearner(Y[(A == 0) & train], X[(A == 0) & train, ],
                        newX = X[test, ], family = 'gaussian', SL.library = sl_lib)
  
  Q0[test] <- Q0.sl$SL.predict
  coefs.Q0 <- coefs.Q0 + Q0.sl$coef

  
  g.sl <- SuperLearner(A[train], X[train, ], newX = X[test, ],
                       family='binomial', SL.library = sl_lib)
  
  g[test] <- g.sl$SL.predict
  coefs.g <- coefs.g + g.sl$coef
}


# Compute ATT Estimate
n1 <- sum(A)
att_hat <- (sum(A * (Y - Q0)) - sum((1 - A) * g / (1 - g) * (Y - Q0))) / n1


# Compute Influence Function and Standard Error

P_A1 <- mean(A)


IF_part1 <- (A / P_A1) * (Y - Q0 - theta_hat)
IF_part2 <- ((1 - A) * g / (P_A1 * (1 - g))) * (Y - Q0)
IF <- IF_part1 - IF_part2

se <- (1/sqrt(N)) * sd(IF)


att_ci_high <- att_hat + 1.96 * se
att_ci_low <- att_hat - 1.96 * se

ci <- theta_hat + c(-1.96, 1.96) * se 

results <- data.frame(
  ATT = theta_hat,
  SE = se,
  CI_low = att_ci_low,
  CI_high = att_ci_high
)
print(results)



# ------------------------------------------------------
# ATT using npcausal
# ------------------------------------------------------

att_npcausal <- att(y = Y, a = A, x = X, nsplits = 5, sl.lib = sl_lib)



# ------------------------------------------------------
# ATT using TMLE
# ------------------------------------------------------



N = length(A)
nsplits = 5
split_inds <- sample(rep(1:nsplits, ceiling(N/2))[1:N])


QA1 <- as_tibble(matrix(0, nrow = N, ncol = 2,
                        dimnames = list(c(), c("Q0", "Q1"))))
Q0 <- rep(NA, N)
g <- rep(NA, N)



coefs.Q0 <- rep(0, length(sl_lib))
coefs.g <- rep(0, length(sl_lib))

# Estimate Nuisance Parameters using Cross-Fitting
for (vfold in 1:nsplits) {
  train <- split_inds != vfold
  test <- split_inds == vfold
  if (nsplits == 1) {
    train <- test
  }
  Q0.sl <- SuperLearner(Y[(A == 0) & train], X[(A == 0) & train, ],
                        newX = X[test, ], family = 'gaussian', SL.library = sl_lib)
  
  Q0[test] <- Q0.sl$SL.predict
  coefs.Q0 <- coefs.Q0 + Q0.sl$coef
  
  
  g.sl <- SuperLearner(A[train], X[train, ], newX = X[test, ],
                       family='binomial', SL.library = sl_lib)
  
  g[test] <- g.sl$SL.predict
  coefs.g <- coefs.g + g.sl$coef
}

# ------------------------------------------------------
# TMLE Targeting Step
# ------------------------------------------------------
P_A1 <- mean(A)
g <- pmin(g, 0.99)
g <- pmax(g, 0.01)

# Construct clever covariate H
H <- ifelse(A == 1, 1 / P_A1, -g / (P_A1 * (1 - g)))

# Fit fluctuation model (update Q0 to Q0_star)
epsilon_model <- glm(Y ~ -1 + H + offset(Q0), family = binomial)
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

# ------------------------------------------------------
# Print Results
# ------------------------------------------------------
results <- data.frame(
  ATT = theta_hat,
  SE = se,
  CI_low = ci[1],
  CI_high = ci[2]
)
print(results)

  