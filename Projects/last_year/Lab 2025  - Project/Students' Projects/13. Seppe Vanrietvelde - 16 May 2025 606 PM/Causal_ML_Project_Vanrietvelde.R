########################################   Question 2	  ######################################
library(ggplot2)
library(hdm)
library(npcausal)
library(SuperLearner)
library(tmle)

set.seed(123)
sl_lib <- c("SL.glm", "SL.ranger", "SL.glmnet")
data(pension)
Y <- pension$net_tfa
A <- pension$e401
X <- pension[, c("age", "inc", "educ", "fsize", "marr", "twoearn", "db", "pira", "hown")]

##########################################   Naive method

model <- lm(Y ~ A, data = pension)
coef(model)["A"] # 19559.34
confint(model)["A",] # (16999.90 22118.79)

##########################################   AIPW - ATE	

aipw_ate <- ate(y = Y, a = A, x = X, nsplits = 5, sl.lib = sl_lib)
saveRDS(aipw_ate, file = "aipw_ate.rds")
#      parameter       est        se     ci.ll     ci.ul pval
# 1      E{Y(0)} 14355.164 1038.4525 12319.797 16390.531    0
# 2      E{Y(1)} 21776.117  833.1272 20143.188 23409.046    0
# 3 E{Y(1)-Y(0)}  7420.953 1243.9635  4982.785  9859.121    0

aipw_ate <- readRDS("aipw_ate.rds")
aipw_ate$res$est[3] # 7420.953
aipw_ate$res[c("ci.ll","ci.ul")][3, ] # 4982.785 9859.121

# Check the propensity score estimates

g1 <- aipw_ate$nuis$pi_a1
ggplot(data.frame(ps = g1, A = factor(A)),
       aes(ps, fill = A, colour = A)) +
  geom_density(alpha = .3) +
  labs(x = "Estimated P(A = 1 | X)", fill = "Treatment")

# Check the weights

w <- A / g1 + (1 - A) / (1 - g1)
summary(w)

hist(w, breaks = 40) # no values with large Weight which might cause unstable results

##########################################   AIPW - ATT

aipw_att <- att(y = Y, a = A, x = X, nsplits = 5, sl.lib = sl_lib)
saveRDS(aipw_att, file = "aipw_att.rds")
#       parameter       est       se    ci.ll    ci.ul pval
# 1      E(Y|A=1) 30347.389 1243.033 27911.05 32783.73    0
# 2   E{Y(0)|A=1} 20762.402 1257.845 18297.03 23227.78    0
# 3 E{Y-Y(0)|A=1}  9584.987 1566.401  6514.84 12655.13    0

aipw_att <- readRDS("aipw_att.rds")
aipw_att$res$est[3] # 9584.987
aipw_att$res[c("ci.ll","ci.ul")][3, ] # 6514.84 12655.13

##########################################   TMLE - ATE	

tmle_ate <- tmle(Y = Y,
                 A = A,
                 W = X,
                 Q.SL.library = sl_lib,
                 g.SL.library = sl_lib,
                 family = 'gaussian',
                 V.Q = 5,
                 V.g = 5)
saveRDS(tmle_ate, file = "tmle_ate.rds")

tmle_ate <- readRDS("tmle_ate.rds")
summary(tmle_ate) # 8117.8 (6187.8 10048)

##########################################   TMLE - ATT

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

# TMLE Targeting Step
P_A1 <- mean(A)
g <- pmin(g, 0.99)
g <- pmax(g, 0.01)

# Construct clever covariate H
H <- ifelse(A == 1, 1 / P_A1, -g / (P_A1 * (1 - g)))

# Fit fluctuation model (update Q0 to Q0_star)
epsilon_model <- glm(Y ~ -1 + H + offset(Q0), family = gaussian)
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

theta_hat # 5505.326 -> I don't think this is correct when comparing with AIPW
ci # (1059.293 9951.359)

########################################   Question 3	  ######################################
library(hdm)
library(SuperLearner)

set.seed(123)
sl_lib <- c("SL.glm", "SL.ranger", "SL.glmnet")
data(pension)
Y <- pension$net_tfa
A <- pension$e401
X <- pension[, c("age", "inc", "educ", "fsize", "marr", "twoearn", "db", "pira", "hown")]

N <- length(A)

g.sl <- SuperLearner(Y = A,
                     X = X,
                     SL.library = sl_lib,
                     family = binomial(),
                     cvControl = list(V=5))
ghat_sl <- g.sl$SL.predict

y.sl <- SuperLearner(Y = Y,
                     X = X,
                     SL.library = sl_lib,
                     family = gaussian(),
                     cvControl = list(V=5))
yhat_sl <- y.sl$SL.predict

cov_aipw_by_hand <- mean((A - ghat_sl)*(Y - yhat_sl))
influence_curve_aipw <- (A - ghat_sl)*(Y - yhat_sl) - cov_aipw_by_hand
cov_aipw_by_hand_se <- (1/sqrt(N)) * sd(influence_curve_aipw)
ci_lower <- cov_aipw_by_hand - 1.96 * cov_aipw_by_hand_se
ci_upper <- cov_aipw_by_hand + 1.96 * cov_aipw_by_hand_se

cov_aipw_by_hand # 1085.769
ci_lower # 737.3289
ci_upper # 1434.208

z <- cov_aipw_by_hand / cov_aipw_by_hand_se
p_value <- 2 * (1 - pnorm(abs(z)))
p_value # 1.01e-09

########################################   Question 4	  ######################################
library(hdm)
library(SuperLearner)
library(rlearner)
library(ggplot2)

set.seed(123)
sl_lib <- list("SL.glm","SL.glmnet","SL.ranger")
data(pension)
Y <- pension$net_tfa
A <- pension$e401
X <- pension[, c("age", "fsize", "inc", "educ", "marr", "twoearn", "db", "ira", "hown")]

# estimating propensity score E[A|L]
sl_g  <- SuperLearner(A, X,
                      family    = binomial(),
                      SL.library= sl_lib,
                      cvControl = list(V = 5)) # 5-fold CV

p_hat <- sl_g$SL.predict  

summary(p_hat)

# Min.   :0.06776  
# 1st Qu.:0.19826  
# Median :0.34877  
# Mean   :0.37176  
# 3rd Qu.:0.51179  
# Max.   :0.95355 

# estimating g-computation E[Y|L]
sl_m  <- SuperLearner(Y = Y, X = X,
                      family    = gaussian(),
                      SL.library= sl_lib,
                      cvControl = list(V = 5))

m_hat <- sl_m$SL.predict 

# rboost
mod_mat <- model.matrix(
  ~ age + inc + educ - 1,     
  data = pension
)

rboost.fit <- rboost(x = mod_mat,
                     w = A,
                     y = Y,
                     k_folds = 5,
                     p_hat = p_hat,
                     m_hat = m_hat)


rboost.est <- predict(rboost.fit, mod_mat)

mean(rboost.est) # 650.7762

rboost.est[0:100]

ggplot(data.frame(tau = rboost.est), aes(tau)) +
  geom_histogram(bins = 30, colour = "white", fill = "steelblue") +
  labs(title = "Histogram of estimated individual effects",
       x = "Estimated CATE", y = "Count") +
  theme_minimal()

ggplot(data.frame(age = mod_mat[,"age"], tau = rboost.est),
       aes(age, tau)) +
  geom_point(alpha = .4) +
  geom_smooth(method = "loess", se = FALSE, colour = "red") +
  labs(title = "Treatment effect as a function of age",
       y = "Estimated CATE", x = "Age") +
  theme_minimal()


ggplot(data.frame(inc = mod_mat[,"inc"], tau = rboost.est),
       aes(inc, tau)) +
  geom_point(alpha = .4) +
  geom_smooth(method = "loess", se = FALSE, colour = "red") +
  labs(title = "Treatment effect as a function of income",
       y = "Estimated CATE", x = "Income") +
  theme_minimal()


ggplot(data.frame(educ = mod_mat[,"educ"], tau = rboost.est),
       aes(educ, tau)) +
  geom_point(alpha = .4) +
  geom_smooth(method = "loess", se = FALSE, colour = "red") +
  labs(title = "Treatment effect as a function of education (yrs)",
       y = "Estimated CATE", x = "Years of education") +
  theme_minimal()

########################################   Question 5	  ######################################
library(ggplot2)
library(hdm)
library(npcausal)
library(SuperLearner)

set.seed(123)
sl_lib <- c("SL.glm", "SL.ranger", "SL.glmnet")
data(pension)

# Prepare training and test sets
K <- 3
group.index <- sample(rep(1:K, each = ceiling(dim(pension)[1]/K)))
group.index <- group.index[1:dim(pension)[1]]

pension.trainA <- pension[group.index == 1,]
pension.trainB <- pension[group.index == 2,]
pension.test <- pension[group.index == 3,]
pension.trainAB <- rbind(pension.trainA, pension.trainB)

x.trainA <- pension.trainA[, c("age", "inc", "educ", "fsize", "marr", "twoearn", "db", "pira", "hown")]  
x.trainB <- pension.trainB[, c("age", "inc", "educ", "fsize", "marr", "twoearn", "db", "pira", "hown")]
x.test <- pension.test[, c("age", "inc", "educ", "fsize", "marr", "twoearn", "db", "pira", "hown")]
x.trainAB <- rbind(x.trainA, x.trainB)

a.trainA <- as.vector(pension.trainA$e401)
a.trainB <- as.vector(pension.trainB$e401)
a.test <- as.vector(pension.test$e401)
a.trainAB <- c(a.trainA, a.trainB)

y.trainA <- as.vector(pension.trainA$net_tfa)
y.trainB <- as.vector(pension.trainB$net_tfa)
y.test <- as.vector(pension.test$net_tfa)
y.trainAB <- c(y.trainA, y.trainB)

# Estimate nuissance parameters
model.ps.trainA <- SuperLearner(Y = a.trainA, 
                                X = x.trainA, 
                                family = binomial(), 
                                cvControl = list(V=5), 
                                SL.library = sl_lib)
model.ps.trainB <- SuperLearner(Y = a.trainB, 
                                X = x.trainB, 
                                family = binomial(), 
                                cvControl = list(V=5), 
                                SL.library = sl_lib)

pension.trainA$ps.SL <- as.vector(predict(model.ps.trainB, x.trainA, onlySL = TRUE)$pred)
pension.trainB$ps.SL <- as.vector(predict(model.ps.trainA, x.trainB, onlySL = TRUE)$pred)

model.Q1.trainA <-  SuperLearner(Y = y.trainA[a.trainA == 1], 
                                 X = x.trainA[a.trainA == 1,],
                                 cvControl = list(V=5),
                                 SL.library = sl_lib)
model.Q1.trainB <-  SuperLearner(Y = y.trainB[a.trainB == 1], 
                                 X = x.trainB[a.trainB == 1,],
                                 cvControl = list(V=5),
                                 SL.library = sl_lib)

pension.trainA$Q1 <- as.vector(predict(model.Q1.trainB, newdata = x.trainA, OnlySL = TRUE)$pred)
pension.trainB$Q1 <- as.vector(predict(model.Q1.trainA, newdata = x.trainB, OnlySL = TRUE)$pred)

# Compute pseudooutcomes
pseudooutcome.DR.trainA <- pension.trainA$Q1 + pension.trainA$e401/pension.trainA$ps.SL * (pension.trainA$net_tfa - pension.trainA$Q1)
pseudooutcome.DR.trainB <- pension.trainB$Q1 + pension.trainB$e401/pension.trainB$ps.SL * (pension.trainB$net_tfa - pension.trainB$Q1)
pseudooutcome.DR.trainAB <- c(pseudooutcome.DR.trainA, pseudooutcome.DR.trainB)

# Estimate counterfactual net financial assets when being eligible for enrolling
model.pseudooutcome.DR <- SuperLearner(Y = pseudooutcome.DR.trainAB, 
                                       X = x.trainAB[, c("age", "inc", "educ")], 
                                       cvControl = list(V=5), 
                                       SL.library = sl_lib)

pension.test$pred.DR <- predict(model.pseudooutcome.DR, newdata = x.test[, c("age", "inc", "educ")], onlySL = TRUE)$pred

# Results
ggplot(data.frame(tau = pension.test$pred.DR), aes(tau)) +
  geom_histogram(bins = 30, colour = "white", fill = "steelblue") +
  labs(x = "Counterfactual net financial assets when being eligible for enrolling", y = "Count") +
  theme_minimal()

# FROM THIS POINT I AM NOT SURE (suggestion of chatgpt)

# Estimate counterfactual prediction error
model.ps.trainAB <- SuperLearner(Y = a.trainAB, 
                                 X = x.trainAB, 
                                 family = binomial(),
                                 cvControl = list(V=5),
                                 SL.library = sl_lib)
pension.test$ps.SL <- as.vector(predict(model.ps.trainAB, newdata = x.test, onlySL = TRUE)$pred)

cf_error <- mean((a.test / pension.test$ps.SL) * (y.test - pension.test$pred.DR)^2)
cf_error # 2778967855

# Compute standard error for counterfactual prediction error using the efficient influence curve
influence_curve <- (a.test / pension.test$ps.SL) * (y.test - pension.test$pred.DR)^2 - cf_error
se <- sd(influence_curve) / sqrt(length(influence_curve))

# Compute the confidence interval for the counterfactual prediction error
ci_lower_IF <- cf_error - 1.96 * se
ci_lower_IF # 1036506632
ci_upper_IF <- cf_error + 1.96 * se
ci_upper_IF # 4521429078
