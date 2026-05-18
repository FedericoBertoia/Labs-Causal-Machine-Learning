#######################################################
######### PROJECT CAUSAL MACHINE LEARNING ############
######################################################

# Sarah De Souter, Charlotte KoekKoek and Eline Vanderpijpen

library(hdm)
library(devtools)
library(npcausal)
library(tmle)
library(ggplot2)
library(rlearner)
library(SuperLearner)

set.seed(1234)

data(pension)
attach(pension)
data <- pension

SL.library <- c("SL.glm", "SL.ranger", "SL.glmnet","SL.gam")

X <- data[, c("age", "inc", "fsize", "educ", "marr", "twoearn", "db", "ira", "hown")]
A <- data[,c("e401")]
Y <- data[,"net_tfa"]
Xr <-  data[,c("age", "inc", "educ")]

####################################################
#####################Question 2#####################
####################################################


##AIPW
#ATE
aipw <- ate(y = Y, a = A, x = X, nsplits = 5, sl.lib = SL.library)
aipw$res
head(aipw$ifvals)

#propensity scores
head(aipw$nuis)
max(aipw$nuis$pi_a1)
min(aipw$nuis$pi_a1)
mean(aipw$nuis$pi_a1)

#plots of propensity scores
a0 <- aipw$nuis$pi_a0
df_a0 <- data.frame(a0)
ggplot(data = df_a0, aes(x=a0)) + geom_histogram(color = "white", fill="pink")
a1 <- aipw$nuis$pi_a1
df_a1 <- data.frame(a1)
ggplot(data = df_a1, aes(x=a1)) + geom_histogram(color = "white", fill="lightblue")

#ATT
aipw2 <- att(y=Y, a=A, x=X, nsplits=5, sl.lib = SL.library)
aipw2$res
head(aipw2$ifvals)

#propensity scores
head(aipw2$nuis)
max(aipw2$nuis$pi)
min(aipw2$nuis$pi)

#plot of propensity scores
a <- aipw2$nuis$pi
df_a <- data.frame(a)
ggplot(data = df_a, aes(x=a)) + geom_histogram(color = "white", fill="pink")

##TMLE
tmle_est <- tmle(Y=Y, A=A, W=X, Q.SL.library = SL.library, g.SL.library = SL.library, V.Q=5,  V.g = 5)
tmle_est
tmle_est$estimates$ATE
tmle_est$estimates$ATT

#propensity scores + plot
propensity_score <- tmle_est$g$g1W
df_propscor <- data.frame(propensity_score)
ggplot(data = df_propscor, aes(x=propensity_score)) + geom_histogram(color = "white", fill="pink")
min(df_propscor)
max(df_propscor)

##No covariate adjustment
#AIPW
#by hand because not supported by package npcausal
N = length(A)
nsplits = 5
split_inds <- sample(rep(1:nsplits, ceiling(N/2))[1:N])

##ATE
g <- rep(NA, N)
Q0 <- rep(NA, N)
Q1 <- rep(NA, N)
bound_pi <- function(pi, lower=0.025, upper=0.975) {
  pi[pi < lower] <- lower
  pi[pi > upper] <- upper
  return(pi)
}

# Estimate Nuisance Parameters using Cross-Fitting
for (vfold in 1:nsplits) {
  train <- split_inds != vfold
  test <- split_inds == vfold
  g[test] <- mean(A[train])
  dat <- as.data.frame(cbind(Y[train], A[train]))
  Qmod <- glm(Y~A, data = dat)
  dat$A <- 0
  Q0[test] <- predict(Qmod, newdata=dat, type = "response")
  dat$A <- 1
  Q1[test] <- predict(Qmod, newdata=dat, type = "response")
}

g1 <- g
g0 <- 1-g
Q <- A*Q1 + (1-A)*Q0

ate_hat <- mean(((A/g1) - (1-A)/g0)*(Y - Q) + Q1 - Q0)

IF <- ((A/g1) - (1-A)/g0)*(Y - Q) + Q1 - Q0 - ate_hat

se <- (1/sqrt(N)) * sd(IF)

ate_ci_high <- ate_hat + 1.96 * se
ate_ci_low <- ate_hat - 1.96 * se

results_ate <- data.frame(
  ATE = ate_hat,
  SE = se,
  CI_low = ate_ci_low,
  CI_high = ate_ci_high
)
results_ate

#ATT
Q0 <- rep(NA, N)
g <- rep(NA, N)


# Estimate Nuisance Parameters using Cross-Fitting
for (vfold in 1:nsplits) {
  train <- split_inds != vfold
  test <- split_inds == vfold
  Q0[test] <- mean(Y[(A == 0) & train])
  g[test] <- mean(A[train])
}

# Compute ATT Estimate
n1 <- sum(A)
att_hat <- (sum(A * (Y - Q0)) - sum((1 - A) * g / (1 - g) * (Y - Q0))) / n1


# Compute Influence Function and Standard Error

P_A1 <- mean(A)
IF <- (A / P_A1) * (Y - Q0 - att_hat) - ((1 - A) * g / (P_A1 * (1 - g))) * (Y - Q0)

se <- (1/sqrt(N)) * sd(IF)


att_ci_high <- att_hat + 1.96 * se
att_ci_low <- att_hat - 1.96 * se

results_att <- data.frame(
  ATT = att_hat,
  SE = se,
  CI_low = att_ci_low,
  CI_high = att_ci_high
)
results_att

#TMLE
N <- matrix(0, nrow= length(Y), ncol = 2)
tmle_est_nocov <- tmle(Y=Y, A=A, W=N, Q.SL.library = SL.library, g.SL.library = SL.library, V.Q=5,  V.g = 5)
tmle_est_nocov
tmle_est_nocov$estimates$ATE
tmle_est_nocov$estimates$ATT

####################################################
#####################Question 3#####################
####################################################

EAL = SuperLearner(A, X, SL.library= SL.library, family='binomial', cvControl = list(5) )
EYL = SuperLearner(Y, X, SL.library = SL.library, family = 'gaussian', cvControl = list(5))
a <- EAL$SL.predict
y <- EAL$SL.predict
theta = mean((A-a)*(Y-y))
theta
N = length(A)
se = (1 / sqrt(N))*sd((A-a)*(Y-y))
se
Z = theta/se
p = 2*(1-pnorm(abs(Z)))
p
lowerbound = theta - 1.96*se
upperbound = theta + 1.96*se
lowerbound
upperbound

####################################################
#####################Question 4#####################
####################################################
set.seed(1234)

# for cross-fitting
splits4 <- sample(rep(1:5, length.out = nrow(data)))

p_hats <- rep(0, length(A))
m_hats <- rep(0, length(A))
for (i in 1:5){
  test_inds <- which(splits4 == i)
  train_inds <- which(splits4 != i)
  
  # models for propensity scores and outcomes using all covariates
  modelp <- SuperLearner(A[train_inds], X[train_inds,], newX = X[test_inds,], family="binomial", SL.library = SL.library, cvControl = list(5))
  p_hat <- modelp$SL.predict 
  p_hats[test_inds] <- p_hat
  
  modelm <- SuperLearner(Y[train_inds], X[train_inds,], newX = X[test_inds,], family="gaussian", SL.library = SL.library, cvControl = list(5))
  m_hat <- modelm$SL.predict
  m_hats[test_inds] <- m_hat
}

# effects in function of age, income and years of education 
rlearner_fit <- rboost(as.matrix(Xr), A, Y, k_folds = 5, p_hat = p_hats, m_hat = m_hats)
rlearner_estimates <- predict(rlearner_fit)

hist_all <- ggplot(as.data.frame(rlearner_estimates), aes(x = rlearner_estimates)) +
  geom_histogram(bins = 20, fill = "grey", color = "black") +
  labs(x = "Effect estimates",
       y = "Count") +
  theme(plot.title = element_text(hjust = 0.5))
hist_all

summary(rlearner_estimates)

# propensity scores:  
eligible <- which(A == 1)
not_eligible <- setdiff(1:length(A), eligible)
groups  <- c(rep("Eligible", length(eligible)),
             rep("Ineligible", length(not_eligible)))
ps <- c(p_hats[eligible], p_hats[not_eligible])
df_ps <- data.frame(Group = groups, Value = ps)
ggplot(df_ps, aes(x = Value, fill = Group)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 30) +
  labs(x = "Propensity Scores",
       y = "Count",
       fill = "Group") +
  scale_fill_manual(values = c("green", "red")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
min(ps)

# rboost using only age
X_age <- as.matrix(X[,"age"])
rlearner_fit_age <- rboost(X_age, A, Y, k_folds = 5, p_hat = p_hats, m_hat = m_hats)
rlearner_estimates_age <- predict(rlearner_fit_age)

hist_age <- ggplot(as.data.frame(rlearner_estimates_age), aes(x = rlearner_estimates_age)) +
  geom_histogram(bins = 20, fill = "lightblue", color = "black") +
  labs(x = "Effect estimates",
       y = "Count") +
  theme(plot.title = element_text(hjust = 0.5))
hist_age

ggplot(as.data.frame(rlearner_estimates_age), aes(x = X_age, y = rlearner_estimates_age)) +
  geom_point(color = "lightblue", size = 2) +
  labs(x = "Age",
       y = "Estimated effect") +
  theme(plot.title = element_text(hjust = 0.5))

summary(rlearner_estimates_age)


# rboost using only income
X_inc <- as.matrix(X[,"inc"])
rlearner_fit_inc <- rboost(X_inc, A, Y, k_folds = 5, p_hat = p_hats, m_hat = m_hats)
rlearner_estimates_inc <- predict(rlearner_fit_inc)

hist_inc <- ggplot(as.data.frame(rlearner_estimates_inc), aes(x = rlearner_estimates_inc)) +
  geom_histogram(bins = 20, fill = "lightgreen", color = "black") +
  labs(x = "Effect estimates",
       y = "Count") +
  theme(plot.title = element_text(hjust = 0.5))
hist_inc

summary(rlearner_estimates_inc)

ggplot(as.data.frame(rlearner_estimates_inc), aes(x = X_inc, y = rlearner_estimates_inc)) +
  geom_point(color = "lightgreen", size = 2) +
  labs(x = "Income",
       y = "Estimated effect") +
  theme(plot.title = element_text(hjust = 0.5))

# can the different effect sizes for higher incomes be explained by one of the confounders? 
# tested with all the confounders but no clear pattern
ggplot(as.data.frame(rlearner_estimates_inc), aes(x = X_inc, y = rlearner_estimates_inc, color = factor(X[,c("fsize")]))) +
  geom_point(size = 2) +
  labs(title = "R learner effect estimates in function of income",
       x = "Income",
       y = "Estimated effect", color = "factor") +
  theme(plot.title = element_text(hjust = 0.5))


# rboost using only education
X_educ <- as.matrix(X[,"educ"])
rlearner_fit_educ <- rboost(X_educ, A, Y, k_folds = 5, p_hat = p_hats, m_hat = m_hats)
rlearner_estimates_educ <- predict(rlearner_fit_educ)

hist_educ <- ggplot(as.data.frame(rlearner_estimates_educ), aes(x = rlearner_estimates_educ)) +
  geom_histogram(bins = 20, fill = "pink", color = "black") +
  labs(x = "Effect estimates",
       y = "Count") +
  theme(plot.title = element_text(hjust = 0.5))
hist_educ

summary(rlearner_estimates_educ)

ggplot(as.data.frame(rlearner_estimates_educ), aes(x = X_educ, y = rlearner_estimates_educ)) +
  geom_point(color = "pink", size = 2) +
  labs(x = "Years of education",
       y = "Estimated effect") +
  theme(plot.title = element_text(hjust = 0.5))

####################################################
#####################Question 5#####################
####################################################

set.seed(1234)

########### Estimation counterfactual net financial assets ################
K <- 5
index1 <- sample(rep(1:K, length.out = nrow(data)))
predictions <- vector(length = nrow(data))
for(i in 1:K){
  test <- which(index1 == i)
  train <- which(index1 != i)
  
  Ytrain <- Y[train]
  Atrain <- A[train]
  Xrtrain <- Xr[train,]
  Xtrain <- X[train,]
  
  index2 <- sample(rep(1:K, length.out = length(train)))
  C <- vector(length = length(train))
  for(j in 1:K){
    test2 <- which(index2 == j)
    train2 <- which(index2 != j)
    Ytrain2 <- Ytrain[train2]
    Atrain2 <- Atrain[train2]
    Xtrain2 <- Xtrain[train2,]
    y1 <- Ytrain2[which(Atrain2 == 1)]
    x1 <- Xtrain2[which(Atrain2 ==1),]
    model.y <- SuperLearner(Y = y1, 
                            X = data.frame(x1), 
                            SL.library = SL.library)
    model.a <- SuperLearner(Y = Atrain2, 
                            X = data.frame(Xtrain2), 
                            SL.library = SL.library)
    newx <- data.frame(Xtrain[test2,])
    pred.ps <- as.vector(predict(model.a,newdata=newx, onlySL = TRUE)$pred)
    pred.y <- as.vector(predict(model.y,newdata=newx, onlySL = TRUE)$pred)
    pseudo.y <- pred.y + Atrain[test2]*(Ytrain[test2]-pred.y)/pred.ps
    C[test2] <- pseudo.y
  }
  model.pseudo <- SuperLearner(Y = C, 
                               X = data.frame(Xrtrain), 
                               SL.library = SL.library)
  predictions[test] <- as.vector(predict(model.pseudo,data.frame(Xr[test,]), onlySL = TRUE)$pred)
}

####### Plots #########
ggplot(as.data.frame(predictions),aes(x=predictions)) +
  geom_histogram(alpha = 0.4, position = "identity", bins = 30, aes(y = after_stat(density)),
                 color = "black")+labs(x = expression(tilde(Y)))

df_y <- data.frame(Y=Y,prediction=predictions)
ggplot(df_y,mapping = aes(x=Y,y=prediction)) + geom_point() + geom_abline(intercept= 0,slope=1,col = 'blue') + labs(x = 'Y',y = expression(tilde(Y)),title = 'All samples')

df_y1 <- data.frame(Y=Y[which(A==1)],prediction=predictions[which(A==1)])
ggplot(df_y1,mapping = aes(x=Y,y=prediction)) + geom_point() + geom_abline(intercept= 0,slope=1,col = 'blue') +labs(x = 'Y',y = expression(tilde(Y)),title = 'Samples in subgroup with E=1')

df_difference <- data.frame(Y=Y,difference=Y-predictions)
ggplot(df_difference,mapping = aes(x=Y,y=difference)) + geom_point()+geom_hline(yintercept = 0,col = 'blue')+labs(x = expression(Y),y = expression(Y-tilde(Y)), title = 'All samples')

df_difference1 <- data.frame(Y=Y[which(A==1)],difference=Y[which(A==1)]-predictions[which(A==1)])
ggplot(df_difference1,mapping = aes(x=Y,y=difference)) + geom_point()+geom_hline(yintercept = 0,col='blue')+labs(x = expression(Y),y = expression(Y-tilde(Y)),title = 'Samples in subgroup with E=1')

########### Estimation counterfactual prediction error #######

Ytilde <- (Y - predictions)^2

aipw3 <- ate(y=Ytilde,a=A,x=data.frame(X),sl.lib= SL.library,nsplits = 5)
aipw3$res

