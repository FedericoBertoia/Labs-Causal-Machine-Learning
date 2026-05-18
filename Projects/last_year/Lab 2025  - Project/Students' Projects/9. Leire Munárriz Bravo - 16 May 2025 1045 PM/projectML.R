#Question 3
N <- dim(pension)[1]
Y <- pension$net_tfa
A <- pension$e401
L <- pension[, c("age", "inc", "fsize", "educ", "marr", 
                 "twoearn", "db", "ira", "hown")]

SL.library <- c( "SL.mean","SL.glm","SL.earth","SL.ranger")


#family = guassian as the variable is continuous
#Y given L
q.model <- SuperLearner(Y = Y, X = L, 
                        SL.library = SL.library, 
                        family = "gaussian", cvControl = list(V = 10))
#A given L
p.model <- SuperLearner(Y = A, X = L, 
                        SL.library = SL.library, 
                        family = "binomial", cvControl = list(V = 10))
p.model

#obtain expectations
q <- q.model$SL.predict
p <- p.model$SL.predict

#compute EIF

theta <- mean((A - p)*(Y - q))
eif<- (A - p)*(Y - q) - theta
se <- sd(eif)/sqrt(length(eif))

ci_lower <- theta - 1.96 * se
ci_upper <- theta + 1.96 * se

#Z-statistic
z_stat <- thetahat / se

#Two-sided p-value
p_value <- 2 * (1 - pnorm(abs(z_stat)))


#Question 5

sl_lib <- c("SL.mean", "SL.glm", "SL.glmnet", "SL.ranger")
nfolds <- 10
split_inds <- generate_folds(length(A), nfolds)

mu_hat <- estimate_mu_nonpar(Y, A, X=L, sl_lib = sl_lib, nsplits = nfolds, family = "gaussian")$preds

#propensity scores per observation
pi_hat <- estimate_pi_nonpar(A, X=L, sl_lib = sl_lib, split_inds = split_inds)$preds
pi_hat <- bound_pi(pi_hat)

C <- mu_hat$mu1 + (Y - mu_hat$mu1)*A/pi_hat

K <- 10
set.seed(123)  
N <- nrow(L)
folds <- sample(rep(1:K, length.out = N))


Y1_tilde <- rep(NA, N)

SL.library <- c("SL.mean", "SL.glm", "SL.earth", "SL.ranger")

#Perform cross-fitting
for (k in 1:K) {
  
  train_idx <- which(folds != k)
  valid_idx <- which(folds == k)

  dr_model <- SuperLearner(Y = C[train_idx],
                           X = L[train_idx, ],
                           newX = L[valid_idx, ],
                           family = gaussian(),
                           SL.library = SL.library)
  
  Y1_tilde[valid_idx] <- dr_model$SL.predict
}

hist(Y1_tilde, breaks = 30, main = "Predicted Counterfactual net_tfa",
     xlab = "Predicted net_tfa under 401(k) Eligibility", col = "lightblue")


squared_errors <- (Y[A == 1] - Y1_tilde[A == 1])^2

# Mean squared error
mse <- mean(squared_errors)

#Standard error
se <- sd(squared_errors) / sqrt(length(squared_errors))

#CI
ci_lower <- mse - 1.96 * se
ci_upper <- mse + 1.96 * se



