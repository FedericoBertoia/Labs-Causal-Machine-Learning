# Libraries ----
library(dplyr)
library(ggplot2)
library(gridExtra)
library(hdm)
library(randomForest)
library(ranger)
library(rlearner) 
library(sandwich)
library(SuperLearner)
library(tibble)
library(tmle)
source("functions.R")


#Q1 ----

# Q2 ----
data("pension")
A <- pension$e401
Y <- pension$net_tfa
X <- pension[,c("age","inc","fsize","marr","educ","twoearn","pira","db", "hown")]

sl_lib <- list("SL.glm","SL.step","SL.glm.interaction", "SL.ranger", "SL.glmnet")

mu_hat_est <- estimate_mu_nonpar(Y = Y,
                                       A = A,
                                       X = X,
                                       sl_lib = sl_lib,
                                       nsplits = 5,
                                       family = 'gaussian')
mu_hat <- mu_hat_est$preds
pi_hat_fit <- estimate_pi_nonpar(A = A,
                                       X = X,
                                       sl_lib = sl_lib,
                                       nsplits = 5)
pi_hat <- pi_hat_fit$preds
pi_hat_bounded <- bound_pi(pi_hat)


ate_aipw_est <- est_aipw(Y, A, pi_hat_bounded, mu_hat)

ate_aipw <- ate_aipw_est[1]; ate_aipw
ate_aipw_se <- ate_aipw_est[2]; ate_aipw_se


ate_tmle_est <- est_tmle(Y, A, X, mu_hat, pi_hat_bounded)

ate_tmle <- ate_tmle_est[1]; ate_tmle
ate_tmle_se <- ate_tmle_est[2]; ate_tmle_se


ate_tmle + 1.96*ate_tmle_se
ate_tmle - 1.96*ate_tmle_se



naive_ate <- mean(Y[A == 1]) - mean(Y[A == 0]); naive_ate
naive_se <- sd(Y[A == 1]) / sqrt(sum(A == 1)) +
  sd(Y[A == 0]) / sqrt(sum(A == 0)) ; naive_se




#Q4 -----
library(doParallel)
registerDoParallel(cores = 4)

r_fit <- rboost(as.matrix(X),A,Y)
saveRDS(r_fit, "r_fit.rds")
tau_hat <- predict(r_fit, as.matrix(X))

hist(tau_hat, main = "Spread of personal treatment effects", xlab = "predicted effect")
p1 <- ggplot(data=NULL, aes(x=X$age, y=tau_hat)) + 
  geom_smooth() +
  labs(title = "Personal treatment effect vs Age", x = "Age", y = "Effect Estimate")
p2 <- ggplot(data=NULL, aes(x=X$inc, y=tau_hat)) + 
  geom_smooth() +
  labs(title = "Personal treatment effect vs income", x = "Income", y = "Effect Estimate")
p3 <- ggplot(data=NULL, aes(x=X$educ, y=tau_hat)) + 
  geom_smooth() +
  labs(title = "Personal treatment effect vs education", x = "Years of education", y = "Effect Estimate")
grid.arrange(p1,p2,p3, nrow=1)
