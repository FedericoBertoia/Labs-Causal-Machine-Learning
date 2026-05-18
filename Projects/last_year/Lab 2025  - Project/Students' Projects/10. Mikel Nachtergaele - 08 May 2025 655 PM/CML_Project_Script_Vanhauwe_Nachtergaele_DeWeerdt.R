###########################Q4
#load libraries
library(hdm)
library(gbm)
library(ggplot2)
library(devtools) 
install_github("xnie/rlearner")
library(rlearner)
library(npcausal)
library(SuperLearner)
library(tmle)
library(haven) 
library(fBasics)

#load data
data(pension)
attach(pension)
print(pension) #Quick check

### Question 2 ====

# Identify required variables
Y <- pension$net_tfa
A <- pension$e401
X <- pension[, c("age", "inc", "fsize", "educ", "marr", 
                 "twoearn", "db", "ira", "hown")]

# Choose libraries of interest for SUperLearner
SL.library <- c("SL.glm", "SL.gam", "SL.glmnet")

# Estimate ATE using AIPW
aipw_ate <- ate(y = Y, a = A, x = X, nsplits = 5, sl.lib = SL.library)
aipw_ate$res

aipw_att <- att(y = Y, a = A, x = X, nsplits = 5, sl.lib = SL.library)
aipw_att$res

# Estimate ATE using TMLE
tmle_results <- tmle(Y = Y, A = A, W = X, 
                     Q.SL.library = SL.library, g.SL.library = SL.library,
                     V.Q = 5, V.g = 5)
summary(tmle_results)

# Estimate unadjusted effect
unadjusted <- lm(net_tfa ~ e401, data = pension)
summary(unadjusted)
confint(unadjusted)

# Ckeck of positivity assumption and extreme cases 

# Estimate propensity scores: P(A = 1 | X)
ps_model <- SuperLearner(Y = A, X = X, 
                         family = binomial(),
                         SL.library = SL.library)

# Extract predicted propensity scores
prop_scores <- ps_model$SL.predict

# Check summary
summary(prop_scores)

# Plot histogram
hist(prop_scores, 
     main = "Histogram of Estimated Propensity Scores", 
     xlab = "Propensity Score", 
     breaks = 50)

# Count extreme values
sum(prop_scores < 0.01)
sum(prop_scores > 0.99)

### Question 3 ====

# Step 1: Estimate E[A | L]
A_fit <- SuperLearner(Y = A, X = X, SL.library = SL.library)
A_pred <- A_fit$SL.predict

# Step 2: Estimate E[Y | L]
Y_fit <- SuperLearner(Y = Y, X = X, SL.library = SL.library)
Y_pred <- Y_fit$SL.predict

# Step 3: Estimate theta
influence_terms <- (A - A_pred) * (Y - Y_pred)
theta_hat <- mean(influence_terms)

# Step 4: Estimate standard error via influence function
n <- nrow(pension)
std_error <- sd(influence_terms) / sqrt(n)

# Step 5: Test H0: theta = 0
z_score <- theta_hat / std_error
p_value <- 2 * (1 - pnorm(abs(z_score)))

# Print results
cat("Estimated theta (E{Cov(A, Y | L)}):", theta_hat, "\n")
cat("Standard Error:", std_error, "\n")
cat("z-score:", z_score, "\n")
cat("p-value:", p_value, "\n")

### Question 4 ====

#outcome
Y <- net_tfa
summary(Y)

#treatment
D <- e401
mean(D)

#model matrix + controls
X <- model.matrix(~ -1 + age + inc + educ + fsize + marr + twoearn + db + pira + 
                    hown, data = pension)

#remove constant variables
X <- X[, which(apply(X, 2, var) != 0)]

#centering
demean <- function(x) x - mean(x)
X <- apply(X, 2, FUN = demean)

#r-learner
rboost_fit <- rboost(X, D, Y, k_folds = 5)

#predict individual treatment effects
rboost_est <- predict(rboost_fit, X)

#add effect estimates back to pension data for analysis
pension$rboost_est <- rboost_est

#histogram of estimated effects
hist(rboost_est, breaks = 50,
     main = "Estimated Effects of 401(k) Eligibility",
     xlab = "Estimated Treatment Effect",
     col = "skyblue", border = "white")

#visualize heterogeneity by Age
ggplot(pension, aes(x = age, y = rboost_est)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", color = "blue", se = FALSE) +
  labs(title = "Estimated Effect vs Age",
       x = "Age", y = "Estimated Effect")

#visualize heterogeneity by Income
ggplot(pension, aes(x = inc, y = rboost_est)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", color = "green", se = FALSE) +
  labs(title = "Estimated Effect vs Income",
       x = "Income", y = "Estimated Effect")

#visualize heterogeneity by Education
ggplot(pension, aes(x = educ, y = rboost_est)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", color = "purple", se = FALSE) +
  labs(title = "Estimated Effect vs Years of Education",
       x = "Years of Education", y = "Estimated Effect")

### Question 5 ====

set.seed(123)

#Relevant columns
columns <- c("e401", "net_tfa", "age", "inc", "fsize", "educ", "marr", "db",
             "pira", "hown")
pension <- pension[, columns]

# Create training and test data sets
# K = k_folds + 1
k_folds <- 5
K <- k_folds + 1
group.index <- sample(rep(1:K, each = ceiling(dim(pension)[1]/K)))
group.index <- group.index[1:dim(pension)[1]]

#pension$group <- group.index

#Control variables
x_col <- c("age", "inc", "fsize", "educ", "marr", "db", "pira", "hown")

x <- pension[, x_col]

a <- as.vector(pension$e401)

y <- as.vector(pension$net_tfa)

#pension$ps.SL[group.index == 1] <- 1
#Superlearner libraries
SL.library <- c("SL.glm", "SL.ranger", "SL.glmnet")

for (i in 1:k_folds){
  # Fit the propensity score model on the training data
  model.ps.traink <- SuperLearner(Y = a[group.index != i], X = x[group.index != i, ],
                                  family = binomial(), SL.library = SL.library)
  # Predict propensity score on left out fold
  pension$ps.SL[group.index == i] <- as.vector(predict(model.ps.traink, x[group.index == i, ], onlySL = TRUE)$pred)
}


for (i in 1:k_folds){
  # Outcome prediction model
  model.Q0.traink <- SuperLearner(Y = y[(group.index != i & a == 0)], X = x[(group.index != i & a==0), ],
                                  SL.library = SL.library)
  model.Q1.traink <- SuperLearner(Y = y[(group.index != i & a == 1)], X = x[(group.index != i & a==1), ],
                                  SL.library = SL.library)
  # Add predicted counterfactual to the training sets. 
  pension$Q0[group.index == i] <- as.vector(predict(model.Q0.traink, x[group.index == i, ], onlySL = TRUE)$pred)
  pension$Q1[group.index == i] <- as.vector(predict(model.Q1.traink, x[group.index == i, ], onlySL = TRUE)$pred)
}

# Compute the pseudooutcomes in the training sets
pension$pseudooutcome.DR <- pension$Q1 + pension$e401/pension$ps.SL * (pension$net_tfa - pension$Q1)

# Regress the pseudooutcomes in the combined training data set
model.pseudooutcome.DR <- SuperLearner(Y = pension$pseudooutcome.DR[group.index != 6], X = x[group.index != 6, ], SL.library = SL.library)

# Predict CATE on the test set
pension$CATE.DR[group.index == 6] <- predict(model.pseudooutcome.DR, newdata = x[group.index == 6, ], onlySL = TRUE)$pred

# Plot the treatment effects obtained on the test set on a histogram.
hist(pension$CATE.DR[group.index == 6], breaks = 50,
     xlab="Net financial assets", 
     main="Estimated counterfactual net financial assets",
      col = "skyblue", border = "white")

median(pension$CATE.DR[group.index == 6])
################################################################################
#Standard errors

N <- dim(pension[group.index == 6, ])[1]

# Fit the propensity score model on the combined training data
model.ps.trainAB <- SuperLearner(Y = a[group.index != 6], X = x[group.index != 6, ], 
                                 family = binomial(), SL.library = SL.library)

# Compute PS on the test set
pension$ps.SL[group.index == 6] <- as.vector(predict(model.ps.trainAB, x[group.index == 6, ], onlySL = TRUE)$pred)

### Outcome model on treated and control
model.Q1.trainAB <- SuperLearner(Y = y[(group.index != 6 & a == 1)], 
                                 X = x[(group.index != 6 & a == 1),], 
                                 SL.library = SL.library)

model.Q0.trainAB <- SuperLearner(Y = y[(group.index != 6 & a == 0)], 
                                 X = x[(group.index != 6 & a == 0),], 
                                 SL.library = SL.library)


#Prediction on test set
pension$Q1[group.index == 6] <- as.vector(predict(model.Q1.trainAB, newdata = x[group.index == 6, ], OnlySL = TRUE)$pred)
pension$Q0[group.index == 6] <- as.vector(predict(model.Q0.trainAB, newdata = x[group.index == 6, ], OnlySL = TRUE)$pred)

###### Compute the pseudooutcomes in the test set
pension$pseudooutcome.DR <- pension$Q1 + pension$e401/pension$ps.SL * (pension$net_tfa - pension$Q1)

y1.aipw <- mean(pension$pseudooutcome.DR[group.index == 6]); y1.aipw

# Compute standard error for ATE using the efficient influence curve
eic.y1.aipw <- pension$pseudooutcome.DR[group.index == 6] - y1.aipw
se.y1.aipw <- (1/sqrt(N)) * sd(eic.y1.aipw); se.y1.aipw

y1_low <- y1.aipw - 1.96 * se.y1.aipw; y1_low
y1_high <- y1.aipw + 1.96 * se.y1.aipw; y1_high

################################################################################

#visualize heterogeneity by Age
ggplot(pension[group.index == 6, ], aes(x = age, y = CATE.DR)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", color = "blue", se = FALSE) +
  labs(title = "Counterfactual net finanial assets vs Age",
       x = "Age", y = "Counterfactual net finanial assets")

#visualize heterogeneity by Income
ggplot(pension[group.index == 6, ], aes(x = inc, y = CATE.DR)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", color = "green", se = FALSE) +
  labs(title = "Counterfactual net finanial assets vs Income",
       x = "Income", y = "Counterfactual net finanial assets")

#visualize heterogeneity by Education
ggplot(pension[group.index == 6, ], aes(x = educ, y = CATE.DR)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", color = "purple", se = FALSE) +
  labs(title = "Counterfactual net finanial assets vs Years of Education",
       x = "Years of Education", y = "Counterfactual net finanial assets")
