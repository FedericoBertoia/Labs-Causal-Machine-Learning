library(hdm)
library(usethis)
library(devtools) 
library(npcausal)
library(SuperLearner)

set.seed(123)
data(pension)
head(pension)

##############
# Question 2 #
##############

table(pension$net_tfa, pension$e401)  

adjusted = net_tfa  ~ e401 + age + inc + fsize + educ + marr + twoearn + db + ira + hown
unadjusted = net_tfa ~ e401

model = lm(unadjusted, data = pension)
summary(model)
confint(model)


Y <- pension$net_tfa
A <- pension$e401
L <- pension[,c("age", "inc", "fsize", "educ", "marr", "twoearn", "db", "ira", "hown")]

library(causaldrf)
sl_lib <- list("SL.glm","SL.randomForest","SL.xgboost")
aipw <- ate(y = Y, a = A, x = L, nsplits = 5, sl.lib = sl_lib) 

aipw$res

library(tmle)

# Create folds for cross-fitting
N <- length(A)
K <- 5
folds <- sample(rep(1:K, length.out = N))

# Initialize vectors to store cross-fitted predictions
Q0 <- Q1 <- g1 <- rep(NA, N)

# Loop over folds
for (k in 1:K) {
  train_idx <- which(folds != k)
  valid_idx <- which(folds == k)
  
  # Fit Q on training set
  SL.Q <- SuperLearner(
    Y = Y[train_idx],
    X = data.frame(A = A[train_idx], L[train_idx, ]),
    SL.library = sl_lib,
    family = gaussian()
  )
  
  # Predict Q1 and Q0 on validation fold
  Q1[valid_idx] <- predict(SL.Q, newdata = data.frame(A = 1, L[valid_idx, ]), onlySL = TRUE)$pred
  Q0[valid_idx] <- predict(SL.Q, newdata = data.frame(A = 0, L[valid_idx, ]), onlySL = TRUE)$pred
  
  # Fit g on training set
  SL.g <- SuperLearner(
    Y = A[train_idx],
    X = L[train_idx, ],
    SL.library = sl_lib,
    family = binomial()
  )
  
  # Predict g1 on validation fold
  g1[valid_idx] <- predict(SL.g, newdata = L[valid_idx, ], onlySL = TRUE)$pred
}

TMLE <- tmle(
  Y = Y,
  A = A,
  W = L,
  Q = cbind(Q0, Q1),
  g1W = g1,
  family = "gaussian"
)
summary(TMLE)


##############
# Question 3 #
##############
# Create folds for cross-fitting
N <- length(A)
K <- 5
folds <- sample(rep(1:K, length.out = N))

# Initialize vectors to store cross-fitted predictions
P <-g1<- rep(NA, N)

# Loop over folds
for (k in 1:K) {
  train_idx <- which(folds != k)
  valid_idx <- which(folds == k)
  
  # Fit P on training set
  SL.P <- SuperLearner(
    Y = Y[train_idx],
    X = L[train_idx, ],
    SL.library = sl_lib,
    family = gaussian()
  )
  
  # Predict P on validation fold
  P[valid_idx] <- predict(SL.P, newdata = L[valid_idx, ], onlySL = TRUE)$pred
  
  # Fit g on training set
  SL.g <- SuperLearner(
    Y = A[train_idx],
    X = L[train_idx, ],
    SL.library = sl_lib,
    family = binomial()
  )
  
  # Predict g1 on validation fold
  g1[valid_idx] <- predict(SL.g, newdata = L[valid_idx, ], onlySL = TRUE)$pred
}

library(dplyr)
#Cov(A,Y∣L)=E[(A−E[A∣L])(Y−E[Y∣L])∣L]
#E{Cov(A,Y∣L)} = E{E[(A−E[A∣L])(Y−E[Y∣L])∣L]}=E[(A−E[A∣L])(Y−E[Y∣L])]

theta_by_hand <- mean((A - g1)*
                      (Y - P)) 


influence_curve_theta <- {A-g1}*{Y-P}
                          - theta_by_hand

theta_hat <- mean(influence_curve_theta)

se_theta <- (1/sqrt(N)) * sd(influence_curve_theta)

theta_hat
se_theta

z_stat <- theta_hat / se_theta
p_value <- 2 * (1 - pnorm(abs(z_stat)))

# 95% CI
ci_ll <- theta_hat - 1.96 * se_theta
ci_ul <- theta_hat + 1.96 * se_theta

ci_ll
ci_ul

z_stat
p_value

# check the stability of inverse probability weights
ipw_weights <- A / g1 + (1 - A) / (1 - g1)
boxplot(ipw_weights, main = "Boxplot of Inverse Probability Weights", ylab = "IPW Weights")

# check  influence functions
boxplot(influence_curve_theta, main = "Boxplot of Influence Function", ylab = "Influence Function")

