library(ggplot2)
library(SuperLearner)
library(npcausal)
library(hdm)
data(pension)

Y <- pension$net_tfa
A <- pension$e401
L <- pension[, c("age", "inc", "fsize", "educ", "marr", "twoearn", "db", "ira", "hown")]

N <- nrow(pension)
######################## Question 2 (Amber) ##############################

SL.library <- c("SL.gam", "SL.randomForest", "SL.glm")

set.seed(2001)

## AIPW
res.ate <- ate(Y, A, L, sl.lib=SL.library, nsplits=5)
res.att <- att(Y, A, L, sl.lib=SL.library, nsplits=5)

res.AIPW <- data.frame(Estimand=c("ATE","ATT"),
                       Estimate=c(res.ate$res[3,2], res.att$res[3,2]),
                       'CI lowerbound'=c(res.ate$res[3,4], res.att$res[3,4]),
                       'CI upperbound'=c(res.ate$res[3,5], res.att$res[3,5]))

## TMLE
# Create 5 folds for cross-fitting
K <- 5
group.index <- sample(rep(1:K, each = ceiling(N/K)))
group.index <- group.index[1:N]
Q0 <- Q1 <- ps <-  rep(NA, N)
for (fold in 1:K){
  #Make initial predictions in each fold, training a model on the rest
  test_id <- which(group.index==fold)
  train_id <- which(group.index!=fold)
  
  y.test <- Y[test_id]
  a.test <- A[test_id]
  l.test <- L[test_id,]
  
  y.train <- Y[train_id]
  a.train <- A[train_id]
  l.train <- L[train_id,]
  
  # Model for outcome
  Q.model <- SuperLearner(Y=y.train,
                          X=data.frame(cbind(A=a.train,l.train)),
                          SL.library=SL.library,
                          family=gaussian())
  # Model for propensity score
  ps.model <- SuperLearner(Y=a.train,
                           X=l.train,
                           SL.library=SL.library,
                           family=binomial())
  
  # Predictions for outcome and propensity scores
  Q0[test_id] <- as.vector(predict(Q.model, newdata=data.frame(A=0,l.test))$pred)
  Q1[test_id] <- as.vector(predict(Q.model, newdata=data.frame(A=1,l.test))$pred)
  ps[test_id] <- as.vector(predict(ps.model, newdata=l.test)$pred)
}
library(tmle)

res.tmle <- tmle(Y=Y,
                 A=A, 
                 W=L,
                 Q=cbind(Q0,Q1),
                 g1W=ps)
summary(res.tmle)

## Unadjusted
unadjusted.ATE <- lm(net_tfa~e401,data=pension)
s <- summary(unadjusted.ATE)
ci <- confint(unadjusted.ATE)
data.frame(Estimand="ATE",
           Estimate=s$coefficients[2,1],
           'CI lowerbound'=ci[2,1],
           'CI upperbound'=ci[2,2])


############################ Question 3 (Feng) #####################################
n <- length(Y)

#5-fold cross-fitting  
k <- 5
est <- rep(NA,n)
folds <- sample(rep(1:k, length.out = n))

for (i in 1:k) {
  train_indices <- which(folds != i)
  pred_indices <- which(folds == i)
  
  A_train <- A[train_indices]
  L_train <- L[train_indices, ]
  Y_train <- Y[train_indices]
  
  SL.library <- c("SL.gam", "SL.randomForest", "SL.glm")
  prop_score_train <- SuperLearner(Y = A_train, X = data.frame(L_train), SL.library = SL.library, family = binomial())
  outcome_reg_train <- SuperLearner(Y = Y_train, X = data.frame(L_train), SL.library = SL.library, family = gaussian())
  
  A_pred <- A[pred_indices]
  L_pred <- L[pred_indices, ]
  Y_pred <- Y[pred_indices]
  
  prop_score_pred <- predict(prop_score_train, newdata = data.frame(L_pred))$pred
  outcome_reg_pred <- predict(outcome_reg_train, newdata = data.frame(L_pred))$pred
  
  est[pred_indices] <-  (A_pred-prop_score_pred)*(Y_pred-outcome_reg_pred)
}

#results
print(mean(est)) #estimation value of theta
print(sd(est)/sqrt(n)) #standard error of theta
print(sqrt(n) * mean(est) / sd(est)) #t test statistic
print(1-pnorm(sqrt(n) * mean(est) / sd(est))) #p-value

############################ Question 4 (Amber) #####################################
Q <- ps <- pseudo_outcomes <- weights <-  rep(NA, N)
for (fold in 1:K){
  # Define training and test/validation folds
  test_id <- which(group.index==fold)
  train_id <- which(group.index!=fold)
  
  y.test <- Y[test_id]
  a.test <- A[test_id]
  l.test <- L[test_id,]
  
  y.train <- Y[train_id]
  a.train <- A[train_id]
  l.train <- L[train_id,]
  
  # Model for outcome
  Q.model <- SuperLearner(Y=y.train,
                          X=l.train,
                          SL.library=SL.library,
                          family=gaussian())
  # Model for propensity score
  ps.model <- SuperLearner(Y=a.train,
                           X=l.train,
                           SL.library=SL.library,
                           family=binomial())
  
  # Predictions for outcome and propensity scores
  Q.fold <- Q[test_id] <-  as.vector(predict(Q.model, l.test)$pred)
  ps.fold <- ps[test_id] <- as.vector(predict(ps.model, newdata=l.test)$pred)
  
  # Compute pseudo_outcomes and weights
  pseudo_outcomes[test_id] <- (y.test-Q.fold)/(a.test-ps.fold)
  weights[test_id] <- (a.test-ps.fold)^2
}

CATE.model <- SuperLearner(Y=pseudo_outcomes,
                           X=L[,c("age","inc","educ")],
                           SL.library=SL.library,
                           family=gaussian(),
                           obsWeights = weights)

CATE.pred <- CATE.model$SL.predict
hist(CATE.pred, nclass=100, xlab="Individual CATE estimates", main="")

R.age.model <- SuperLearner(Y=pseudo_outcomes,
                            X=data.frame(L[,c("age")]), 
                            SL.library=SL.library,
                            family=gaussian())
R.age <- R.age.model$SL.predict

R.income.model <- SuperLearner(Y=pseudo_outcomes,
                               X=data.frame(L[,c("inc")]), 
                               SL.library=SL.library,
                               family=gaussian())
R.income <- R.income.model$SL.predict

R.years.model <- SuperLearner(Y=pseudo_outcomes,
                              X=data.frame(L[,c("educ")]), 
                              SL.library=SL.library,
                              family=gaussian())
R.years <- R.years.model$SL.predict

predictions <- data.frame(cbind(pension, pred.age=R.age, pred.income=R.income, pred.years=R.years))


ggplot(predictions, aes(x=age, y=pred.age))+geom_point()+
  labs(x="Age",y="Individual effect") +geom_smooth(method="lm", se=F)
ggplot(predictions, aes(x=inc, y=pred.income))+geom_point()+
  labs(x="Income",y="Individual effect") +geom_smooth(method="lm",se=F)
ggplot(predictions, aes(x=educ, y=pred.years))+geom_point()+
  labs(x="Years of education",y="Individual effect") +geom_smooth(method="lm",se=F)
################################# Question 5 #####################################

L_select <- cbind(pension$age, pension$inc, pension$educ)

#3-fold cross-fitting  
#step 1：
k <- 3
folds <- sample(rep(1:k, length.out = n))
pseudo_outcome <- rep(NA, n)
CATE_pred <- rep(NA, n)

#step 2:
for (i in 1:k) {
  
  last_indices <- which(folds == i)
  L_select_last <- L_select[last_indices, ]
  
  for (j in (1:k)[-i]) {
    train_indices <- which(folds == j)
    pred_indices <- which(folds==(1:k)[-c(i,j)])
    
    A_train <- A[train_indices]
    L_train <- L[train_indices, ]
    Y_train <- Y[train_indices]
    
    L_train_treat <- L_train[A_train==1, ]
    Y_train_treat <- Y_train[A_train==1]
    
    prop_score_train <- SuperLearner(Y = A_train, X = data.frame(L_train), SL.library = SL.library, family = binomial())
    outcome_reg_train_treat <- SuperLearner(Y = Y_train_treat, X = data.frame(L_train_treat), SL.library = SL.library, family = gaussian())
    
    A_pred <- A[pred_indices]
    L_pred <- L[pred_indices, ]
    Y_pred <- Y[pred_indices]
    
    prop_score_pred <- predict(prop_score_train, newdata = data.frame(L_pred))$pred
    outcome_reg_pred_treat <- predict(outcome_reg_train_treat, newdata = data.frame(L_pred))$pred
    
    pseudo_outcome[pred_indices] <-  outcome_reg_pred_treat+A_pred*(Y_pred-outcome_reg_pred_treat)/prop_score_pred
  }
  #pseudo outcome C_i
  pseudo_outcome <- pseudo_outcome[c(train_indices, pred_indices)]
  L_select_comb <- L_select[c(train_indices, pred_indices), ]
  
  #step 3:
  CATE_train <- SuperLearner(Y = pseudo_outcome, X = data.frame(L_select_comb), SL.library = SL.library, family = gaussian())
  
  # tilde{Y}^1
  CATE_pred[last_indices] <- predict(CATE_train, newdata = data.frame(L_select_last))$pred
  
}
hist(CATE_pred, main = "Histogram of DR-learner", xlab = "DR-learner", ylab = "Frequency")     

#estimate counterfactual prediction error
Z <- (Y-CATE_pred)^2

#5-fold cross-fitting  
k <- 5
pred_error <- rep(NA,n)
folds <- sample(rep(1:k, length.out = n))

for (i in 1:k) {
  train_indices <- which(folds != i)
  pred_indices <- which(folds == i)
  
  A_train <- A[train_indices]
  L_train <- L[train_indices, ]
  Z_train <- Z[train_indices]
  
  L_train_treat <- L_train[A_train==1, ]
  Z_train_treat <- Z_train[A_train==1]
  
  prop_score_train <- SuperLearner(Y = A_train, X = data.frame(L_train), SL.library = SL.library, family = binomial())
  outcome_reg_train_treat <- SuperLearner(Y = Z_train_treat, X = data.frame(L_train_treat), SL.library = SL.library, family = gaussian())
  
  A_pred <- A[pred_indices]
  L_pred <- L[pred_indices, ]
  Z_pred <- Z[pred_indices]
  
  prop_score_pred <- predict(prop_score_train, newdata = data.frame(L_pred))$pred
  outcome_reg_pred_treat <- predict(outcome_reg_train_treat, newdata = data.frame(L_pred))$pred
  
  pred_error[pred_indices] <-  outcome_reg_pred_treat+A_pred*(Z_pred-outcome_reg_pred_treat)/prop_score_pred
}
#results
print(mean(pred_error))#estimation value
print(mean(pred_error)-qnorm(0.975)*sd(pred_error)/sqrt(n))#lower bound of confidence interval
print(mean(pred_error)+qnorm(0.975)*sd(pred_error)/sqrt(n))#upper bound of confidence interval