###########################################################################
######################          PROJECT CML          ######################         
###########################################################################

library(hdm)
library(SuperLearner)
library(twang)
library(tibble)
library(dplyr)
data(pension)

setwd("C:/Users/milvervo/OneDrive - UGent/Documents/UGent/causal machine learning/project")
#setwd("C:/Users/mavergou/OneDrive - UGent/Documents/Courses/Causal Machine Learning/Project")

source("functions.R")

###########################################################################
######################          QUESTION 1           ######################         
###########################################################################

# See report

###########################################################################
######################          QUESTION 2           ######################         
###########################################################################

set.seed(9000)


###################################
#####       Functions         #####
###################################

# Function to identify outlier rows
remove_outlier_rows <- function(vec, mat) {
  # Sanity checks
  if (length(vec) != nrow(mat)) {
    stop("Length of vector must match number of rows in the matrix.")
  }
  
  # Compute IQR and outlier bounds
  Q1 <- quantile(vec, 0.25, na.rm = TRUE)
  Q3 <- quantile(vec, 0.75, na.rm = TRUE)
  IQR <- Q3 - Q1
  lower_bound <- Q1 - 1.5 * IQR
  upper_bound <- Q3 + 1.5 * IQR
  
  # Logical vector for non-outliers
  non_outlier_mask <- vec >= lower_bound & vec <= upper_bound
  
  # Return filtered matrix
  return(mat[non_outlier_mask, , drop = FALSE])
}


###########################################################################
######################            AIPW               ######################         
###########################################################################

Y <- pension$net_tfa            # Net financial assets
A <- pension$e401               # Eligibility for 401(k)
X <- pension[, c("age", "inc", "fsize", "educ", "marr", "twoearn", "db", "pira", "hown")]


###############################
#####         ATT         #####
###############################

library(devtools) #install.packages("devtools")
library(npcausal) #install_github("ehkennedy/npcausal")

sl_lib <- list("SL.glm", "SL.ranger")

## 1. ATT (Main result)
# Calculate ATT: libary produces honest predictions, as predictions are made on a holdout sample not used for training. 
# Source code here: https://github.com/ehkennedy/npcausal/blob/master/R/att.R
att <- att(y = Y, a = A, x = X, nsplits = 5, sl.lib = sl_lib); att$res
pi_hat = att$nuis$pi 

# Visualize propensity scores
Apihat = as.data.frame(cbind(A, pi_hat))

pi_A1 <- Apihat$pi_hat[Apihat$A == 1]
pi_A0 <- Apihat$pi_hat[Apihat$A == 0]

tail(sort(pi_A1))
tail(sort(pi_A0))

head(sort(pi_A0))
head(sort(pi_A1))

plot(density(pi_A1, na.rm = TRUE), col = "red", lwd = 2, 
     main = "Kernel Density Estimation of propensity score", 
     xlab = "Nonparametric propensity score", ylim = c(0, max(density(pi_A1)$y, density(pi_A0)$y)))
lines(density(pi_A0, na.rm = TRUE), col = "blue", lwd = 2)
legend("topright", legend = c("A = 1", "A = 0"), col = c("red", "blue"), lwd = 2)

# Summarize influence functions
head(att$ifvals)
summary(att$ifvals)

## 2. ATT (No outliers)
# Robustness: redo ATT with excluding outliers of influence curve 
influence_scores_ate <- rowSums(att$ifvals)

filtered_X <- remove_outlier_rows(att$ifvals$V2, X)
filtered_A <- remove_outlier_rows(att$ifvals$V2, as.matrix(A))
filtered_Y <- remove_outlier_rows(att$ifvals$V2, as.matrix(Y))

filtered_att <- att(y = filtered_Y, a = filtered_A, x = filtered_X, nsplits = 5, sl.lib = sl_lib); filtered_att$res

X_naive <- data.frame(intercept = rep(1, length(Y)))

## 3. ATT (Naive, no covariate adjustment - This section is based on efficient influence curves, not the package)
# Naive (using EICs)
n <- length(Y)
p <- mean(A)
EY0 <- mean(Y[A == 0])  # E[Y | A = 0]
EY1 <- mean(Y[A == 1])  # E[Y | A = 1]
EY <- mean(Y)  # E[Y]

naive_att = mean(A/p*(Y-EY0) - (1-A)/(1-p) * p/p * (Y - EY0))/mean(A/p); naive_att
eic = A/p*(Y-EY0 - naive_att) - (1-A)/(1-p) * p/p * (Y - EY0)
naive_se = 1/sqrt(n) * sd(eic); naive_se

# Confidence interval
ATT.AIPW.naive.CI.low <- naive_att - 1.96 * naive_se; ATT.AIPW.naive.CI.low
ATT.AIPW.naive.CI.high <- naive_att + 1.96 * naive_se; ATT.AIPW.naive.CI.high
z_score = naive_att * sqrt(N) * (1/sd(eic)); z_score
p_value = 2 * (1 - pnorm(abs(z_score))); p_value
###############################
#####         ATE         #####
###############################

## 1. ATE (Main result)
# Calculate the ATE
ate <- ate(y = Y, a = A, x = X, nsplits = 5, sl.lib = sl_lib); ate$res

# Summarize influence functions
head(ate$ifvals)

## 2. ATE (No outliers)
# Redo ATE without outliers influence curve 
influence_scores_ate <- rowSums(ate$ifvals)

filtered_X_ate <- remove_outlier_rows(influence_scores_ate, X)
filtered_A_ate <- remove_outlier_rows(influence_scores_ate, as.matrix(A))
filtered_Y_ate <- remove_outlier_rows(influence_scores_ate, as.matrix(Y))

filtered_ate <- ate(y = filtered_Y_ate, a = filtered_A_ate, x = filtered_X_ate, nsplits = 5, sl.lib = sl_lib); filtered_ate$res

## 3. ATE (Naive, no covariate adjustment - This section is based on efficient influence curves, not the package)
# Naive  (using EICs)
naive_ate = mean((A/p - (1-A)/(1-p))*(Y-EY) + EY1 - EY0); naive_ate
eic = (A/p - (1-A)/(1-p))*(Y-EY) + EY1 - EY0 - naive_ate
naive_se = 1/sqrt(n) * sd(eic); naive_se

# Confidence interval
ATE.AIPW.naive.CI.low <- naive_ate - 1.96 * naive_se; ATE.AIPW.naive.CI.low
ATE.AIPW.naive.CI.high <- naive_ate + 1.96 * naive_se; ATE.AIPW.naive.CI.high
z_score = naive_ate * sqrt(N) * (1/sd(eic)); z_score
p_value = 2 * (1 - pnorm(abs(z_score))); p_value

###########################################################################
######################            TMLE               ######################         
###########################################################################

set.seed(9000)


###############################
#####         ATE         #####
###############################

sl_lib <- list("SL.glm", "SL.ranger")

## 1. ATE (Main result)
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
                        newX = X[test, ], family = 'gaussian', SL.library = sl_lib) # Modified for continuous outcome
  
  Q[test, "Q0"] <- Q0.sl$SL.predict
  coefs.Q <- coefs.Q + Q0.sl$coef
  
  Q1.sl <- SuperLearner(Y[(A == 1) & train], X[(A == 1) & train, ],
                        newX = X[test, ], family = 'gaussian', SL.library = sl_lib) # Modified for continuous outcome
  
  Q[test, "Q1"] <- Q1.sl$SL.predict
  coefs.Q <- coefs.Q + Q1.sl$coef
  
  g.sl <- SuperLearner(A[train], X[train, ], newX = X[test, ],
                       family='binomial', SL.library = sl_lib)
  
  g[test] <- g.sl$SL.predict
  coefs.g <- coefs.g + g.sl$coef
}

coefs.Q <- coefs.Q/(2*nsplits); coefs.Q # SuperLearner weights (not used)
coefs.g <- coefs.g/nsplits; coefs.g # SuperLearner weights (not used)

H1 <- A/g
H0 <- (1-A)/(1-g)

# Estimate parameters delta in the fluctuation model and target the predictions
Q.tmle <- as_tibble(matrix(0, nrow = N, ncol = 2,
                           dimnames = list(c(), c("Q0.tmle", "Q1.tmle"))))


delta1.tmle <- coef(glm(Y ~ -1 + offset(Q$Q1) + H1, family = 'gaussian'))
Q.tmle["Q1.tmle"]<- Q$Q1 + delta1.tmle/g

delta0.tmle <- coef(glm(Y ~ -1 + offset(Q$Q0) + H0, family = 'gaussian'))
Q.tmle["Q0.tmle"]<- Q$Q0 + delta0.tmle/(1-g)


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

## 2. ATE (No covariate adjustment)
## TMLE NAIVE
g = mean(A)
Q0 <- mean(Y[A == 0])
Q1 <- mean(Y[A == 1])
Q0 <- rep(Q0, length(Y))
Q1 <- rep(Q1, length(Y))
H1 <- A/g
H0 <- (1-A)/(1-g)

# Estimate parameters delta in the fluctuation model and target the predictions
Q.tmle <- as_tibble(matrix(0, nrow = N, ncol = 2,
                           dimnames = list(c(), c("Q0.tmle", "Q1.tmle"))))


delta1.tmle <- coef(glm(Y ~ -1 + offset(Q1) + H1, family = 'gaussian'))
Q.tmle["Q1.tmle"]<- Q1 + delta1.tmle/g

delta0.tmle <- coef(glm(Y ~ -1 + offset(Q0) + H0, family = 'gaussian'))
Q.tmle["Q0.tmle"]<- Q0 + delta0.tmle/(1-g)

EY1.tmle <- mean(Q.tmle$Q1.tmle); EY1.tmle
EY0.tmle <- mean(Q.tmle$Q0.tmle); EY0.tmle
ATE.tmle <- EY1.tmle - EY0.tmle; ATE.tmle

ATE.tmle <- EY1.tmle - EY0.tmle; ATE.tmle

se.ATE.tmle <- (1/sqrt(N)) * sd(EIC.EY1.tmle-EIC.EY0.tmle); se.ATE.tmle

ATE.tmle.CI.low <- ATE.tmle - 1.96 * se.ATE.tmle; print(ATE.tmle.CI.low)
ATE.tmle.CI.high <- ATE.tmle + 1.96 * se.ATE.tmle; print(ATE.tmle.CI.high)

###############################
#####         ATT         #####
###############################
set.seed(9000)

## 1. ATT (Main Result)
N = length(A)
nsplits = 5
split_inds <- sample(rep(1:nsplits, ceiling(N/2))[1:N])

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


P_A1 <- mean(A)
g <- pmin(g, 0.99)
g <- pmax(g, 0.01)

# Construct clever covariate H
H <- ifelse(A == 1, 1 / P_A1, -g / (P_A1 * (1 - g)))

# Fit fluctuation model (update Q0 to Q0_star)
epsilon_model <- glm(Y ~ -1 + H + offset(Q0), family = "gaussian")
epsilon <- coef(epsilon_model)["H"]
Q0_tmle <- Q0 + epsilon * H

# Compute ATT
ATT.tmle <- mean( (Y[A == 1] - Q0_tmle[A == 1]) ); ATT.tmle

# Influence function and standard error
IF_part1 <- (A / P_A1) * (Y - Q0_tmle - ATT.tmle)
IF_part2 <- ((1 - A) * g / (P_A1 * (1 - g))) * (Y - Q0_tmle)
IF <- IF_part1 - IF_part2
se.ATT.tmle <- (1/sqrt(N)) * sd(IF); se.ATT.tmle

# Confidence interval
ATT.tmle.CI.low <- ATT.tmle - 1.96 * se.ATT.tmle; ATT.tmle.CI.low
ATT.tmle.CI.high <- ATT.tmle + 1.96 * se.ATT.tmle; ATT.tmle.CI.high
z_score = ATT.tmle * sqrt(N) * (1/sd(IF)); z_score
p_value = 2 * (1 - pnorm(abs(z_score))); p_value

## 2. ATT (No covariate adjustment)

## TMLE NAIVE
g = mean(A)
Q0 <- mean(Y[A == 0])
Q1 <- mean(Y[A == 1])
Q0 <- rep(Q0, length(Y))
Q1 <- rep(Q1, length(Y))

H <- ifelse(A == 1, 1 / g, -g / (g * (1 - g)))

# Fit fluctuation model (update Q0 to Q0_star)
epsilon_model <- glm(Y ~ -1 + H + offset(Q0), family = "gaussian")
epsilon <- coef(epsilon_model)["H"]
Q0_tmle <- Q0 + epsilon * H

# Compute ATT
ATT.tmle <- mean( (Y[A == 1] - Q0_tmle[A == 1]) ); ATT.tmle

# Influence function and standard error
IF_part1 <- (A / g) * (Y - Q0_tmle - ATT.tmle)
IF_part2 <- ((1 - A) * g / (g * (1 - g))) * (Y - Q0_tmle)
IF <- IF_part1 - IF_part2
se.ATT.tmle <- (1/sqrt(N)) * sd(IF); se.ATT.tmle

# Confidence interval
ATT.tmle.CI.low <- ATT.tmle - 1.96 * se.ATT.tmle; ATT.tmle.CI.low
ATT.tmle.CI.high <- ATT.tmle + 1.96 * se.ATT.tmle; ATT.tmle.CI.high
z_score = ATT.tmle * sqrt(N) * (1/sd(IF)); z_score
p_value = 2 * (1 - pnorm(abs(z_score))); p_value


###########################################################################
######################          QUESTION 3           ######################         
###########################################################################

set.seed(9000)

###### 1. Outcome regression (this estimates nonparametrically E(Y|L))
y_L <- CV.SuperLearner(Y = Y, X = X, SL.library = sl_lib, V = 5)
y_L_hat <- y_L$SL.predict # these predictions are made within on a holdout sample not used for training

##### 2. Estimator and standard error
N <- length(A)
est_nonpar <- mean((A - pi_hat)*(Y - y_L_hat)); est_nonpar
est_nonpar_eic <- (A - pi_hat)*(Y - y_L_hat) - est_nonpar
est_nonpar_se <- (1/sqrt(N)) * sd(est_nonpar_eic); est_nonpar_se

z_score = est_nonpar * sqrt(N) * (1/sd(est_nonpar_eic)); z_score
ul = est_nonpar + 1.96 * est_nonpar_se; ul
ll = est_nonpar - 1.96 * est_nonpar_se; ll
p_value = 2 * (1 - pnorm(abs(z_score))); p_value


###################################################################################
######################          QUESTION 4                   ######################         
###################################################################################

set.seed(2000)
cate_results <- data.frame()

K <- 5  # 5-fold sample-splitting
group.index <- sample(rep(1:K, length.out = dim(pension)[1]))

# Prepare a data frame to store CATE predictions from both folds
cate_results <- data.frame()

# Loop over the splits
for (i in 1:5) {
  
  # Split data into training and testing based on the current fold
  data.train <- pension[group.index != i,] 
  data.test <- pension[group.index == i,]
  
  # Extract features for both training and testing
  x.train <- data.train[, c("age", "inc", "fsize", "educ", "marr", "twoearn", "db", "pira", "hown")]
  x.test <- data.test[, c("age", "inc", "fsize", "educ", "marr", "twoearn", "db", "pira", "hown")]
  
  a.train <- as.vector(data.train$e401)
  a.test <- as.vector(data.test$e401)
  
  y.train <- as.vector(data.train$net_tfa)
  y.test <- as.vector(data.test$net_tfa)
  
  ##############################
  ####  Nuisance parameters ####
  ##############################
  
  # Outcome prediction model using SuperLearner
  model.Q.train <- CV.SuperLearner(Y = y.train, X = x.train, SL.library = sl_lib, V = 5)  # Continuous outcome
  
  # Compute predictions from the outcome model for the training set
  data.train$Q <- model.Q.train$SL.predict
  
  # Propensity score model for the training set
  model.ps.train <- CV.SuperLearner(Y = a.train, X = x.train, SL.library = sl_lib, family = "binomial", V = 5) # Binary outcome
  
  # Compute predictions from the propensity score model
  data.train$ps.SL <- model.ps.train$SL.predict
  
  #######################################
  #### Pseudooutcomes and prediction ####
  #######################################
  
  # Compute the pseudooutcomes in the training set
  pseudooutcome.R.train <- (data.train$net_tfa - data.train$Q) / (data.train$e401 - data.train$ps.SL)
  
  # Define weights for regression
  w.train <- (a.train - c(data.train$ps.SL))^2
  
  # Use a subset of features for regression
  x.train_short <- x.train[, c("age", "inc", "educ")]
  
  # Fit model using the pseudooutcomes
  model.pseudooutcome.R <- SuperLearner(Y = pseudooutcome.R.train, X = x.train_short, SL.library = c("SL.glm"), obsWeights = w.train, cvControl = list(V = 5))
  
  # Predict CATE on the test set using the trained model
  x.test_short <- x.test[, c("age", "inc", "educ")]
  data.test$CATE.R <- predict(model.pseudooutcome.R, newdata = as.matrix(x.test_short), onlySL = FALSE)$pred
  
  # Store the predictions and corresponding covariates
  fold_results <- data.frame(fold = i, age = x.test$age, inc = x.test$inc, educ = x.test$educ, CATE.R = data.test$CATE.R)
  cate_results <- rbind(cate_results, fold_results)
}

#########################
####  Visualization  ####
#########################


hist(cate_results$CATE.R, xlab="R-learner", main="", breaks = 30)

summary(cate_results$CATE.R)

# Plot all CATE predictions together (from both folds)
library(ggplot2)

ggplot(cate_results, aes(age, CATE.R)) +
  geom_point() + 
  geom_smooth(method = "loess") + 
  ggtitle("Age effect") + 
  xlab("Age") + 
  ylab("CATE")

ggplot(cate_results, aes(inc, CATE.R)) +
  geom_point() + 
  geom_smooth(method = "loess") + 
  ggtitle("Income effect") + 
  xlab("Income") + 
  ylab("CATE")

ggplot(cate_results, aes(educ, CATE.R)) +
  geom_point() + 
  geom_smooth(method = "loess") + 
  ggtitle("Schooling effect") + 
  xlab("Years of education") + 
  ylab("CATE")



###################################################################################
######################          QUESTION 5                   ######################         
###################################################################################

set.seed(2000)
cate_results <- data.frame()
cate_results_test <- data.frame()

K <- 2  # 2-fold sample-splitting
group.index <- sample(rep(1:K, each = ceiling(dim(pension)[1] / K)))
group.index <- group.index[1:dim(pension)[1]]

# Initialize library for SuperLearner
sl_lib <- list("SL.glm", "SL.ranger")

# Loop over the splits
for (i in 1:K) {
  
  ####################### PART 1   #######################
  #####          Prediction of Y_1 hat               ##### 
  
  # Split data into training and testing based on the current fold
  data.train <- pension[group.index == i,] 
  data.test <- pension[group.index != i,]
  
  # Extract features for both training and testing
  x.train <- data.train[, c("age", "inc", "fsize", "educ", "marr", "twoearn", "db", "pira", "hown")]  
  x.test <- data.test[, c("age", "inc", "fsize", "educ", "marr", "twoearn", "db", "pira", "hown")] 
  
  a.train <- as.vector(data.train$e401)
  a.test <- as.vector(data.test$e401)
  
  y.train <- as.vector(data.train$net_tfa)
  y.test <- as.vector(data.test$net_tfa)
  
  # Estimate nuisance parameters (outcome model for A = 1)
  N <- length(a.train)
  nsplits <- 5
  split_inds <- sample(rep(1:nsplits, ceiling(N/2))[1:N])
  
  Q <- as_tibble(matrix(0, nrow = N, ncol = 1, dimnames = list(c(), c("Q1"))))
  coefs.Q <- rep(0, length(sl_lib))
  
  for (vfold in 1:nsplits) {
    train <- split_inds != vfold
    test <- split_inds == vfold
    if (nsplits == 1) {
      train <- test
    }
    
    Q1.sl <- SuperLearner(y.train[(a.train == 1) & train], 
                          x.train[(a.train == 1) & train, ], 
                          newX = x.train[test, ], 
                          family = 'gaussian', 
                          SL.library = sl_lib)  # Modified for continuous outcome
    
    Q[test, "Q1"] <- Q1.sl$SL.predict
    coefs.Q <- coefs.Q + Q1.sl$coef
  }
  
  # Store outcome predictions for the training set
  data.train$Q1 <- Q$Q1
  
  # Propensity score model for the training set
  model.ps.train <- CV.SuperLearner(Y = a.train, X = x.train, SL.library = sl_lib, family = "binomial", V = 5)
  
  # Compute predictions from the propensity score model
  data.train$ps.SL <- model.ps.train$SL.predict
  
  # Compute pseudooutcomes
  pseudooutcome.DR.train <- data.train$Q1 + a.train / data.train$ps.SL * (y.train - data.train$Q)
  
  # Regress the pseudooutcomes in the combined training data set
  x.train_short <- x.train[, c("age", "inc", "educ")]
  model.pseudooutcome.DR <- SuperLearner(Y = pseudooutcome.DR.train, X = x.train_short, SL.library = "SL.glm")
  
  # Predict Y_1hat on the test set
  x.test_short <- x.test[, c("age", "inc", "educ")]
  data.test$CF1 <- predict(model.pseudooutcome.DR, newdata = x.test_short, onlySL = TRUE)$pred
  
  # Store the results for this fold
  folddr_results <- data.frame(fold = i, age = x.test$age, inc = x.test$inc, educ = x.test$educ, CF1.DR = data.test$CF1)
  cate_results <- rbind(cate_results, folddr_results)
  
  
  ####################### PART 2 #######################
  ###### Prediction of pseudo-outcomes for test ########
  N = length(a.test)
  nsplits = 5
  split_inds <- sample(rep(1:nsplits, ceiling(N/2))[1:N])
  
  Q <- as_tibble(matrix(0, nrow = N, ncol = 1,
                        dimnames = list(c(), c("Q1"))))
  
  coefs.Q <- rep(0, length(sl_lib))
  
  for (vfold in 1:nsplits) {
    train <- split_inds != vfold
    test <- split_inds == vfold
    if (nsplits == 1) {
      train <- test
    }
    Q1.sl <- SuperLearner(y.test[(a.test == 1) & train], x.test[(a.test == 1) & train, ],
                          newX = x.test[test, ], family = 'gaussian', SL.library = sl_lib) # Modified for continuous outcome
    
    Q[test, "Q1"] <- Q1.sl$SL.predict
    coefs.Q <- coefs.Q + Q1.sl$coef
    
  }
  
  data.test$Q1 <- Q$Q1
  
  # Propensity score model (only on test set)
  model.ps.test <- CV.SuperLearner(Y = a.test, X = x.test, SL.library = sl_lib, family = "binomial", V=5)
  
  # Compute predictions from the propensity score model
  data.test$ps.SL <- model.ps.test$SL.predict
  
  # Compute pseudooutcomes
  pseudooutcome.DR.test <- data.test$Q1 + a.test/data.test$ps.SL * (y.test - data.test$Q1)
  folddr_results_test <- data.frame(fold = i, age = x.test$age, inc = x.test$inc, educ = x.test$educ, pseudo.DR = pseudooutcome.DR.test)
  cate_results_test <- rbind(cate_results_test, folddr_results_test)
  
}

summary(cate_results$CF1.DR)

# For part 1: Plot the treatment effects obtained on the test set on a histogram
hist(cate_results$CF1.DR, xlab="DR-learner (Y_1 tilde)", main="")

# For part 2: Compute the counterfactual prediction errors
hist(cate_results_test$pseudo.DR, xlab="DR-learner (Y_1)", main="")

filter_idx <- cate_results_test$pseudo.DR >= -50000 & cate_results_test$pseudo.DR <= 250000
hist(cate_results_test$pseudo.DR[filter_idx], xlab="DR-learner (Y_1)", main="")

pred_error <- mean((cate_results$CF1 - cate_results_test$pseudo.DR)^2) ; pred_error
sqrt(pred_error)

pred_error <- mean((cate_results$CF1[filter_idx] - cate_results_test$pseudo.DR[filter_idx])^2) ; pred_error
sqrt(pred_error)
