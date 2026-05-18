### 0 Setup

library(hdm)           # pension data
library(SuperLearner)  # ensemble nuisance estimation
library(tmle)          # TMLE implementation
library(ggplot2)       # plotting
library(sandwich)      # robust SEs

remotes::install_github("ehkennedy/npcausal")
library(npcausal)

#Load some functions from the Lab 5 solutions:

# Propensity-score cross-fitting 
estimate_pi_nonpar <- function(A, X, sl_lib, nsplits = 5, family = binomial()) {
  N <- length(A)
  folds <- sample(rep(seq_len(nsplits), length.out = N))
  preds <- numeric(N)
  coef_sum <- rep(0, length(sl_lib))
  for (v in seq_len(nsplits)) {
    train <- which(folds != v)
    test  <- which(folds == v)
    fit <- SuperLearner(
      Y           = A[train],
      X           = X[train, , drop=FALSE],
      newX        = X[test,  , drop=FALSE],
      SL.library  = sl_lib,
      family      = family,
      method      = "method.NNloglik"
    )
    preds[test]    <- fit$SL.predict
    coef_sum       <- coef_sum + fit$coef
  }
  list(preds = preds, coefs = coef_sum / nsplits)
}

estimate_mu_nonpar <- function(Y, A, X, sl_lib, nsplits = 5, family = gaussian()) {
  N <- length(Y)
  folds <- sample(rep(seq_len(nsplits), length.out = N))
  mu0 <- numeric(N); mu1 <- numeric(N)
  coef0 <- coef1 <- rep(0, length(sl_lib))
  for (v in seq_len(nsplits)) {
    train <- which(folds != v)
    test  <- which(folds == v)
    # Fit on treated = 1
    fit1 <- SuperLearner(
      Y           = Y[A[train]==1 & train],
      X           = X[A[train]==1 & train, , drop=FALSE],
      newX        = X[test, , drop=FALSE],
      SL.library  = sl_lib,
      family      = family
    )
    mu1[test] <- fit1$SL.predict
    coef1     <- coef1 + fit1$coef
    # Fit on treated = 0
    fit0 <- SuperLearner(
      Y           = Y[A[train]==0 & train],
      X           = X[A[train]==0 & train, , drop=FALSE],
      newX        = X[test, , drop=FALSE],
      SL.library  = sl_lib,
      family      = family
    )
    mu0[test] <- fit0$SL.predict
    coef0     <- coef0 + fit0$coef
  }
  list(mu0  = mu0,
       mu1  = mu1,
       coefs = list(coef0 / nsplits, coef1 / nsplits))
}

### 1 Load data and define objects

data(pension)

Y      <- pension$net_tfa
A      <- pension$e401
L      <- pension[, c("age","inc","fsize","educ","marr","twoearn","db","ira","hown")]
sl_lib <- c("SL.glmnet", "SL.gam")

### 2 Unadjusted ATE 

unadj_ATE <- mean(Y[A==1]) - mean(Y[A==0])
unadj_SE  <- sqrt(var(Y[A==1])/sum(A==1) + var(Y[A==0])/sum(A==0))
cat("Unadjusted ATE =", round(unadj_ATE,1), "SE =", round(unadj_SE,1), "\n")

### 3. AIPW via npcausal::ate()

aipw_res <- ate(
  y       = Y,
  a       = A,
  x       = L,
  nsplits = 5,
  sl.lib  = sl_lib,
  trim    = NULL
)
print(aipw_res$res)  # E[Y(0)], E[Y(1)], ATE , SE & 95% CI
cat("\n")

### 4 Manual cross-fitting + AIPW 

# 4a) Propensity scores g(L)
ghat <- estimate_pi_nonpar(A, L, sl_lib, nsplits=5)$preds

# 4b) Outcome regressions Q0(L), Q1(L)
mu_fit <- estimate_mu_nonpar(Y, A, L, sl_lib, nsplits=5)
Q0     <- mu_fit$mu0
Q1     <- mu_fit$mu1
tau    <- Q1 - Q0

# 4c) Influence-curve and ATE
N      <- length(Y)
IC     <- (A/ghat - (1-A)/(1-ghat)) * (Y - ifelse(A==1, Q1, Q0)) + tau
ate_m  <- mean(IC)
se_m   <- sd(IC)/sqrt(N)
cat("Manual AIPW  ATE =", round(ate_m,1), "SE =", round(se_m,1), "\n\n")

### 5 Diagnostics for AIPW

# 5a. Density of g?(L) by A
df_propensity_scores <- data.frame(
  ghat = ghat,
  A    = factor(A, levels = c(0,1), labels = c("Ineligible","Eligible"))
)
ggplot(df_propensity_scores, aes(x = ghat, color = A)) +
  geom_density() +
  labs(
    title = "Density of Propensity Scores g?(L) by Eligibility",
    x     = expression(hat(g)(L)),
    color = "Eligibility"
  )

# 5b. Boxplot of inverse-propensity weights
df_weights <- data.frame(
  w = ifelse(A == 1, 1/ghat, 1/(1 - ghat)),
  A = factor(A, levels = c(0,1), labels = c("Ineligible","Eligible"))
)
ggplot(df_weights, aes(x = A, y = w)) +
  geom_boxplot() +
  labs(
    title = "Inverse-Propensity Weights by Eligibility",
    x     = "Eligibility",
    y     = expression(1 / hat(g)(L) ~ "or" ~ 1/(1 - hat(g)(L)))
  )

summary(IC)

### 6 Covariates and weights

# 1. Build DF with covariates, propensity & weight
df <- data.frame(
  L,
  ghat   = ghat,
  weight = ifelse(A == 1, 1/ghat, 1/(1 - ghat))
)

# 2. Restrict to the eligible (A=1) households
df_e <- subset(df, A == 1)

# 3a. See the top 5 highest?weight treated units
top5 <- df_e[order(-df_e$weight), ][1:5, ]
print(top5)

# 3b. Quick correlations between each L_j and the weight
cors <- sapply(names(L), function(var){
  cor(df_e[[var]], df_e$weight)
})
sort(cors, decreasing=TRUE)

### 7 AIPW ATT via npcausal::att() 

att_res <- att(
  y       = Y,
  a       = A,
  x       = L,
  nsplits = 5,
  sl.lib  = sl_lib)

### 8 TMLE ATE

tmle_fit <- tmle(
  Y               = Y,
  A               = A,
  W               = L,
  Q.SL.library    = sl_lib,
  g.SL.library    = sl_lib,
  family          = "gaussian",
  V.Q             = 5,
  V.g             = 5
)

# Extract and print ATE
ate_tmle    <- tmle_fit$estimates$ATE$psi
se_ate_tmle <- sqrt(tmle_fit$estimates$ATE$var.psi)

### 9 TMLE ATT

# 1) Obtain fitted g(L) and Q0(L)
ghat  <- estimate_pi_nonpar(A, L, sl_lib, nsplits = 5)$preds
mu0_1 <- estimate_mu_nonpar(Y, A, L, sl_lib, nsplits = 5)$mu0

# 2) Fluctuation model on the control group
fluct_mod <- lm(
  Y[A == 0] ~ -1 +
    offset(mu0_1[A == 0]) +
    I(ghat[A == 0] / (1 - ghat[A == 0]))
)
epsilon <- coef(fluct_mod)

# 3) Update Q1 via targeted fluctuation
Q1_star <- mu0_1 + epsilon * (ghat / (1 - ghat))

# 4) Compute ATT and its SE via the influence curve
pA       <- mean(A)
ATT_tmle <- mean((A / pA) * (Y - Q1_star))
IC_att   <- (A / pA) * (Y - Q1_star - ATT_tmle) -
  ((1 - A) * ghat / (1 - ghat)) / pA * (Y - Q1_star)
se_att_tmle <- sd(IC_att) / sqrt(length(Y))

cat(sprintf("TMLE ATT = %.1f (SE = %.1f)\n",
            ATT_tmle, se_att_tmle))

### 9 E{Cov(A,Y|L)}=0

# 7a) Crossfit E[Y | L] with SuperLearner

#Copy-paste lab 5 but without A 
estimate_m_nonpar <- function(Y, X, sl_lib, nsplits = 5, family = gaussian()) {
  N <- length(Y)
  folds <- sample(rep(seq_len(nsplits), length.out = N))
  preds <- numeric(N)
  coef_sum <- rep(0, length(sl_lib))
  for (v in seq_len(nsplits)) {
    train <- which(folds != v)
    test  <- which(folds == v)
    fit <- SuperLearner(
      Y          = Y[train],
      X          = X[train, , drop = FALSE],
      newX       = X[test,  , drop = FALSE],
      SL.library = sl_lib,
      family     = family
    )
    preds[test]  <- fit$SL.predict
    coef_sum     <- coef_sum + fit$coef
  }
  list(preds = preds, coefs = coef_sum / nsplits)
}

# (re)compute g(L) if needed
g_out <- estimate_pi_nonpar(A, L, sl_lib, nsplits = 5)
ghat  <- g_out$preds

# compute m(L) = E[Y | L]
m_out <- estimate_m_nonpar(Y, L, sl_lib, nsplits = 5)
mhat  <- m_out$preds

# form Di = (A - g(L)) * (Y - m(L))
D     <- (A - ghat) * (Y - mhat)
N     <- length(D)

# estimate theta and its SE
theta_hat <- mean(D)
se_hat    <- sd(D) / sqrt(N)

# two-sided p-value for H0: theta = 0
z_stat <- theta_hat / se_hat
pval   <- 2 * (1 - pnorm(abs(z_stat)))