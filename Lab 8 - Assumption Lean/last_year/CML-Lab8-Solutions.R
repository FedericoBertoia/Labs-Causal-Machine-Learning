library(lmtest)
library(sandwich)
library(SuperLearner)

# Load data
file_path <- 
df <- read.table(file_path)
colnames(df) <- c("gender", "plural", "age", "race", "parity", "married", "bwt", "smokeN", "drinkN", "firstep", "welfare", "smoker", "drinker", "wpre", "wgain", "education", "gestation")

############################################################################

# King County 2001 Birth Data
# =========================================================================
#   
# This data file contains outcomes on 2500 children that were born in
# King County in 2001.  
# 
# The variables are (in column order):
#   
# "gender"        M = male, F = female baby
# "plural"        1 = singleton, 2 = twin, 3 = triplet
# "age"           mother's age in years
# "race"          race categories (for mother)
# "parity"        number of previous live born infants
# "married"       Y = yes, N = no
# "bwt"           birth weight in grams
# "smokeN"        number of cigarettes smoked per day during pregnancy
# "drinkN"        number of alcoholic drinks per week during pregnancy
# "firstep"       1 = participant in program; 0 = did not participate
# "welfare"       1 = participant in public assistance program; 0 = did not 
# "smoker"        Y = yes, N = no, U = unknown
# "drinker"       Y = yes, N = no, U = unknown 
# "wpre"          mother's weight in pounds prior to pregnancy
# "wgain"         mother's weight gain in pounds during pregnancy
# "education"     highest grade completed (add 12 + 1 / year of college)
# "gestation"     weeks from last menses to birth of child

############################################################################

##########################
# Setup
##########################

set.seed(123)

# Number of study participants
N <- dim(df)[1]

# We select the subset of relevant variables (in particular we wish to exclude post-treatment variables)
data <- subset(df, select = c(bwt, firstep, gender, age, race, parity, married, smoker, drinker, wpre, education))

# Dichotomize the outcome variable
data$low_bwt <- ifelse(data$bwt < 2500, 1, 0)

# Set the SuperLearner library
SL.library <- c("SL.glm", "SL.ranger")



##########################
# Question 1
##########################

Y <- data$low_bwt
A <- data$firstep
L <- subset(data, select = -c(low_bwt, bwt, firstep))

q.model <- SuperLearner(Y = Y, X = L, SL.library = SL.library, family = "binomial", cvControl = list(V = 5))
p.model <- SuperLearner(Y = A, X = L, SL.library = SL.library, family = "binomial", cvControl = list(V = 5))

q <- q.model$SL.predict
p <- p.model$SL.predict

# --- LM with sandwich SE ---
model <- lm(I(Y-q) ~ -1 + I(A-p))
lm_result <- coeftest(model, vcov = sandwich)

beta_lm <- lm_result[1, 1]
se_lm <- lm_result[1, 2]

beta_lm_lb <- beta_lm - 1.96 * se_lm
beta_lm_ub <- beta_lm + 1.96 * se_lm
  
# --- EIF ---
beta_eif <- mean((A - p) * (Y - q)) / mean((A - p)^2)
se_eif <- (1 / sqrt(N)) * sd((A - p) * (Y - q - beta_eif * (A - p)) / mean((A - p)^2))

beta_eif_lb <- beta_eif - 1.96 * se_eif
beta_eif_ub <- beta_eif + 1.96 * se_eif

# --- AIPW ---
library(npcausal) 
ate.aipw <- ate(y = Y, a = A, x = L, nsplits = 5, sl.lib = SL.library)

ate_aipw <- ate.aipw$res$est[3]
se_aipw <- ate.aipw$res$se[3]

ate_aipw_lb <- ate.aipw$res$ci.ll[3]
ate_aipw_ub <- ate.aipw$res$ci.ul[3]

# --- TMLE ---
library(tmle)
ate.tmle <- tmle(Y = Y, A = A, W = L, Q.SL.library = SL.library, g.SL.library = SL.library, family = 'binomial', V.Q = 5, V.g = 5)

ate_tmle <- ate.tmle$estimates$ATE$psi
se_tmle <- sqrt(ate.tmle$estimates$ATE$var.psi)

ate_tmle_lb <- ate.tmle$estimates$ATE$CI[1]
ate_tmle_ub <- ate.tmle$estimates$ATE$CI[2]



# --- Combine results into a table ---
results <- tibble::tibble(
  Method = c("Linear Model", "Eff. Influence Function", "AIPW", "TMLE"),
  Estimate = c(beta_lm, beta_eif, ate_aipw, ate_tmle),
  StdError = c(se_lm, se_eif, se_aipw, se_tmle),
  CI_Lower = c(beta_lm_lb, beta_eif_lb, ate_aipw_lb, ate_tmle_lb),
  CI_Upper = c(beta_lm_ub, beta_eif_ub, ate_aipw_ub, ate_tmle_ub)
)

print(results)




##########################
# Question 2
##########################
age <- L$age

age_summary <- mean((A - p)^2 * age) / mean((A - p)^2)
se_age_summary <- (1 / sqrt(N)) * sd( (A - p)^2 * (age - age_summary) / mean(p*(1-p)) )

# age_summary2 <- mean((A*(1-p)^2 + (1-A)*p^2) * age) / mean((A*(1-p)^2 + (1-A)*p^2))
# se_age_summary2 <- (1 / sqrt(N)) * sd( (A*(1-p)^2 + (1-A)*p^2) * (age - age_summary) / mean(p*(1-p)) ) 

# age_summary3 <- mean(((p*(1-p) + (1-(2*p))*(A-p)) * age)) / mean(p*(1-p) + (1-(2*p))*(A-p))
# se_age_summary3 <- (1 / sqrt(N)) * sd( ((p*(1-p) + (1-(2*p))*(A-p))) * (age - age_summary) / mean(p*(1-p)) ) 

L$gender <- ifelse(L$gender == 'M', 1, 0)
L$smoker <- ifelse(L$smoker == 'Y', 1, 0)
L$drinker <- ifelse(data$drinker == 'Y', 1, 0)


# Weighted Summary

weighted_summary <- matrix(NA, nrow = length(L), ncol = 4)
rownames(weighted_summary) <- names(L)
colnames(weighted_summary) <- c("weighted_mean", "weighted_se", "weighted_CI_lb", "weighted_CI_ub")

for (varname in names(L)) {
  X <- L[[varname]]
  
  if (is.numeric(X)) {  
    weights <- (A - p)^2
    
    # Weighted mean, SE and CIs
    weighted_mean <- mean(weights * X) / mean(weights)
    
    eif <- weights * (X - weighted_mean) / mean(p*(1-p))
    weighted_se <- (1 / sqrt(N)) * sd(eif)
    
    weighted_CI_lb <- weighted_mean - 1.96 * weighted_se
    weighted_CI_ub <- weighted_mean + 1.96 * weighted_se
    
    # Save results
    weighted_summary[varname, "weighted_mean"] <- weighted_mean
    weighted_summary[varname, "weighted_se"] <- weighted_se
    weighted_summary[varname, "weighted_CI_lb"] <- weighted_CI_lb
    weighted_summary[varname, "weighted_CI_ub"] <- weighted_CI_ub
    
  } else {
    warning(paste("Variable", varname, "is not numeric. Skipping."))
  }
}

weighted_summary <- round(weighted_summary, 3)

print(weighted_summary)


# Classic Summary

classic_summary <- matrix(NA, nrow = length(L), ncol = 4)
rownames(classic_summary) <- names(L)
colnames(classic_summary) <- c("classic_mean", "classic_se", "classic_CI_lb", "classic_CI_ub")

for (varname in names(L)) {
  X <- L[[varname]]
  
  if (is.numeric(X)) {  
    # Classic (unweighted) mean, SE and CIs
    classic_mean <- mean(X)
    classic_se <- sd(X) / sqrt(N)
    
    classic_CI_lb <- classic_mean - 1.96*classic_se
    classic_CI_ub <- classic_mean + 1.96*classic_se
    
    # Save results
    classic_summary[varname, "classic_mean"] <- classic_mean
    classic_summary[varname, "classic_se"] <- classic_se
    classic_summary[varname, "classic_CI_lb"] <- classic_CI_lb
    classic_summary[varname, "classic_CI_ub"] <- classic_CI_ub
    
  } else {
    warning(paste("Variable", varname, "is not numeric. Skipping."))
  }
}

classic_summary <- round(classic_summary, 3)

print(classic_summary)


##########################
# Question 3
##########################

# For the purpose of the analysis we exclude information about parity and participation in the FirstStep programme as adjusting for these variables could remove part of the age effect. 
y <- data$low_bwt
a <- data$age
l <- subset(data, select = -c(low_bwt, bwt, firstep, age, parity))

q.model2 <- SuperLearner(Y = y, X = l, SL.library = SL.library, family = "binomial", cvControl = list(V = 5))
p.model2 <- SuperLearner(Y = a, X = l, SL.library = SL.library, family = "gaussian", cvControl = list(V = 5))

q2 <- q.model2$SL.predict
p2 <- p.model2$SL.predict


# --- LM with sandwich SE ---
model2 <- lm(I(y-q2) ~ -1 + I(a-p2))
lm_result2 <- coeftest(model2, vcov = sandwich)

beta_lm2 <- lm_result2[1, 1]
se_lm2 <- lm_result2[1, 2]

beta_lm_lb2 <- beta_lm2 - 1.96 * se_lm2
beta_lm_ub2 <- beta_lm2 + 1.96 * se_lm2

# --- EIF ---
beta_eif2 <- mean((a - p2) * (y - q2)) / mean((a - p2)^2)
se_eif2 <- (1 / sqrt(N)) * sd((a - p2) * (y - q2 - beta_eif2 * (a - p2)) / mean((a - p2)^2))

beta_eif_lb2 <- beta_eif2 - 1.96 * se_eif2
beta_eif_ub2 <- beta_eif2 + 1.96 * se_eif2


# --- Combine results into a table ---
results2 <- tibble::tibble(
  Method = c("Linear Model", "Eff. Influence Function"),
  Estimate = c(beta_lm2, beta_eif2),
  StdError = c(se_lm2, se_eif2),
  CI_Lower = c(beta_lm_lb2, beta_eif_lb2),
  CI_Upper = c(beta_lm_ub2, beta_eif_ub2)
)

print(results2)


##########################
# Question 4 
##########################

SL.library2 <- c("SL.glm", "SL.randomForest", "SL.earth")

# Split data into 2 parts
I <-  sort(sample(1:N, N/2))
IC <- setdiff(1:N,I)

Y <- data$low_bwt
A <- data$firstep
L <- subset(data, select = -c(low_bwt, bwt, firstep))

# --- FirstSteps on Birth Weight ---

q3I.model <- SuperLearner(Y = Y[I], X = L[I,], SL.library = SL.library2, family = "binomial", cvControl = list(V = 5))
q3IC.model <- SuperLearner(Y = Y[IC], X = L[IC,], SL.library = SL.library2, family = "binomial", cvControl = list(V = 5))

p3I.model <- SuperLearner(Y = A[I], X = L[I,], SL.library = SL.library2, family = "binomial", cvControl = list(V = 5))
p3IC.model <- SuperLearner(Y = A[IC], X = L[IC,], SL.library = SL.library2, family = "binomial", cvControl = list(V = 5))

q3I <- predict(q3IC.model, newdata = L[I,], onlySL = TRUE)$pred
q3IC <- predict(q3I.model, newdata = L[IC,], onlySL = TRUE)$pred

p3I <- predict(p3IC.model, newdata = L[I,], onlySL = TRUE)$pred
p3IC <- predict(p3I.model, newdata = L[IC,], onlySL = TRUE)$pred

q3 <- c(q3I, q3IC)
p3 <- c(p3I, p3IC)

# --- LM with sandwich SE ---
model3 <- lm(I(Y-q3) ~ -1 + I(A-p3))
lm_result3 <- coeftest(model3, vcov = sandwich)

beta_lm3 <- lm_result3[1, 1]
se_lm3 <- lm_result3[1, 2]

beta_lm_lb3 <- beta_lm3 - 1.96 * se_lm3
beta_lm_ub3 <- beta_lm3 + 1.96 * se_lm3

# --- EIF ---
beta_eif3 <- mean((A - p3) * (Y - q3)) / mean((A - p3)^2)
se_eif3 <- (1 / sqrt(N)) * sd((A - p3) * (Y - q3 - beta_eif3 * (A - p3)) / mean((A - p3)^2))

beta_eif_lb3 <- beta_eif3 - 1.96 * se_eif3
beta_eif_ub3 <- beta_eif3 + 1.96 * se_eif3

# --- Combine results into a table ---
results3 <- tibble::tibble(
  Method = c("Linear Model", "Eff. Influence Function"),
  Estimate = c(beta_lm3, beta_eif3),
  StdError = c(se_lm3, se_eif3),
  CI_Lower = c(beta_lm_lb3, beta_eif_lb3),
  CI_Upper = c(beta_lm_ub3, beta_eif_ub3)
)

print(results3)

# --- Age on Birth Weight ---

q4I.model <- SuperLearner(Y = y[I], X = l[I,], SL.library = SL.library2, family = "binomial", cvControl = list(V = 5))
q4IC.model <- SuperLearner(Y = y[IC], X = l[IC,], SL.library = SL.library2, family = "binomial", cvControl = list(V = 5))

p4I.model <- SuperLearner(Y = a[I], X = l[I,], SL.library = SL.library2, family = "gaussian", cvControl = list(V = 5))
p4IC.model <- SuperLearner(Y = a[IC], X = l[IC,], SL.library = SL.library2, family = "gaussian", cvControl = list(V = 5))

q4I <- predict(q4IC.model, newdata = l[I,], onlySL = TRUE)$pred
q4IC <- predict(q4I.model, newdata = l[IC,], onlySL = TRUE)$pred

p4I <- predict(p4IC.model, newdata = l[I,], onlySL = TRUE)$pred
p4IC <- predict(p4I.model, newdata = l[IC,], onlySL = TRUE)$pred

q4 <- c(q4I, q4IC)
p4 <- c(p4I, p4IC)

# --- LM with sandwich SE ---
model4 <- lm(I(Y-q4) ~ -1 + I(A-p4))
lm_result4 <- coeftest(model4, vcov = sandwich)

beta_lm4 <- lm_result4[1, 1]
se_lm4 <- lm_result4[1, 2]

beta_lm_lb4 <- beta_lm4 - 1.96 * se_lm4
beta_lm_ub4 <- beta_lm4 + 1.96 * se_lm4

# --- EIF ---
beta_eif4 <- mean((A - p4) * (Y - q4)) / mean((A - p4)^2)
se_eif4 <- (1 / sqrt(N)) * sd((A - p4) * (Y - q4 - beta_eif4 * (A - p4)) / mean((A - p4)^2))

beta_eif_lb4 <- beta_eif4 - 1.96 * se_eif4
beta_eif_ub4 <- beta_eif4 + 1.96 * se_eif4

# --- Combine results into a table ---
results4 <- tibble::tibble(
  Method = c("Linear Model", "Eff. Influence Function"),
  Estimate = c(beta_lm4, beta_eif4),
  StdError = c(se_lm4, se_eif4),
  CI_Lower = c(beta_lm_lb4, beta_eif_lb4),
  CI_Upper = c(beta_lm_ub4, beta_eif_ub4)
)

print(results4)


##########################
# Questions 5 and 6 (First Steps Program)
##########################

# --- Preliminaries ---

AL <- subset(data, select = -c(low_bwt, bwt))

Y_hat.model <- SuperLearner(Y = Y, X = AL, SL.library = SL.library, family = "binomial", cvControl = list(V = 5))
Y_hat <- Y_hat.model$SL.predict

#################
# Risk Ratio
#################

log_Y <- log(Y_hat)

q_log.model <- SuperLearner(Y = log_Y, X = L, SL.library = SL.library, family = "gaussian", cvControl = list(V = 5))
q_log <- q_log.model$SL.predict

# --- LM with sandwich SE ---
pseudo_Y_log <- log_Y - q_log + ((Y - Y_hat)/(Y_hat))

log_model <- lm(I(pseudo_Y_log) ~ -1 + I(A-p))
log_lm_result <- coeftest(log_model, vcov = sandwich)

log_beta_lm <- log_lm_result[1, 1]
log_se_lm <- log_lm_result[1, 2]

log_beta_lm_lb <- log_beta_lm - 1.96 * log_se_lm
log_beta_lm_ub <- log_beta_lm + 1.96 * log_se_lm

RR_coef <- exp(log_beta_lm)
RR_se <- log_se_lm

RR_CI_lb <- exp(log_beta_lm_lb)
RR_CI_ub <- exp(log_beta_lm_ub)

# --- Combine results into a table ---
RR_results <- tibble::tibble(
  Method = c("Risk Ratio"),
  Estimate = c(RR_coef),
  StdError = c(RR_se),
  CI_Lower = c(RR_CI_lb),
  CI_Upper = c(RR_CI_ub)
)

print(RR_results)

##################
# Odds Ratio
##################

logit_Y <- log(Y_hat / (1 - Y_hat))

q_logit.model <- SuperLearner(Y = logit_Y, X = L, SL.library = SL.library, family = "gaussian", cvControl = list(V = 5))
q_logit <- q_logit.model$SL.predict

# --- LM with sandwich SE ---
pseudo_Y_logistic <- logit_Y - q_logit + ((Y - Y_hat)/(Y_hat*(1-Y_hat)))

logistic_model <- lm(I(pseudo_Y_logistic) ~ -1 + I(A-p))
logistic_lm_result <- coeftest(logistic_model, vcov = sandwich)

logistic_beta_lm <- logistic_lm_result[1, 1]
logistic_se_lm <- logistic_lm_result[1, 2]

logistic_beta_lm_lb <- logistic_beta_lm - 1.96 * logistic_se_lm
logistic_beta_lm_ub <- logistic_beta_lm + 1.96 * logistic_se_lm

OR_coef <- exp(logistic_beta_lm)
OR_se <- logistic_se_lm

OR_CI_lb <- exp(logistic_beta_lm_lb)
OR_CI_ub <- exp(logistic_beta_lm_ub)
  
# --- Combine results into a table ---
OR_results <- tibble::tibble(
  Method = c("Odds Ratio"),
  Estimate = c(OR_coef),
  StdError = c(OR_se),
  CI_Lower = c(OR_CI_lb),
  CI_Upper = c(OR_CI_ub)
)

print(OR_results)


##########################
# Questions 5 and 6 (Age)
##########################


# --- Preliminaries ---

al <- subset(data, select = -c(low_bwt, bwt, firstep, age, parity))

y_hat.model <- SuperLearner(Y = y, X = al, SL.library = SL.library, family = "binomial", cvControl = list(V = 5))
y_hat <- y_hat.model$SL.predict

#################
# Risk Ratio
#################

log_y <- log(y_hat)

q2_log.model <- SuperLearner(Y = log_y, X = l, SL.library = SL.library, family = "gaussian", cvControl = list(V = 5))
q2_log <- q2_log.model$SL.predict

# --- LM with sandwich SE ---
pseudo_y_log <- log_y - q2_log + ((y - y_hat)/(y_hat))

log_model2 <- lm(I(pseudo_y_log) ~ -1 + I(A-p2))
log_lm_result2 <- coeftest(log_model2, vcov = sandwich)

log_beta_lm2 <- log_lm_result2[1, 1]
log_se_lm2 <- log_lm_result2[1, 2]

log_beta_lm_lb2 <- log_beta_lm2 - 1.96 * log_se_lm2
log_beta_lm_ub2 <- log_beta_lm2 + 1.96 * log_se_lm2

RR_coef2 <- exp(log_beta_lm2)
RR_se2 <- log_se_lm2

RR_CI_lb2 <- exp(log_beta_lm_lb2)
RR_CI_ub2 <- exp(log_beta_lm_ub2)

# --- Combine results into a table ---
RR_results2 <- tibble::tibble(
  Method = c("Risk Ratio Age"),
  Estimate = c(RR_coef2),
  StdError = c(RR_se2),
  CI_Lower = c(RR_CI_lb2),
  CI_Upper = c(RR_CI_ub2)
)

print(RR_results2)

##################
# Odds Ratio
##################

logit_y <- log(y_hat / (1 - y_hat))

q2_logit.model <- SuperLearner(Y = logit_y, X = l, SL.library = SL.library, family = "gaussian", cvControl = list(V = 5))
q2_logit <- q2_logit.model$SL.predict

# --- LM with sandwich SE ---
pseudo_y_logistic <- logit_y - q2_logit + ((y - y_hat)/(y_hat*(1-y_hat)))

logistic_model2 <- lm(I(pseudo_y_logistic) ~ -1 + I(A-p))
logistic_lm_result2 <- coeftest(logistic_model2, vcov = sandwich)

logistic_beta_lm2 <- logistic_lm_result2[1, 1]
logistic_se_lm2 <- logistic_lm_result2[1, 2]

logistic_beta_lm_lb2 <- logistic_beta_lm2 - 1.96 * logistic_se_lm2
logistic_beta_lm_ub2 <- logistic_beta_lm2 + 1.96 * logistic_se_lm2

OR_coef2 <- exp(logistic_beta_lm2)
OR_se2 <- logistic_se_lm2

OR_CI_lb2 <- exp(logistic_beta_lm_lb2)
OR_CI_ub2 <- exp(logistic_beta_lm_ub2)

# --- Combine results into a table ---
OR_results2 <- tibble::tibble(
  Method = c("Odds Ratio Age"),
  Estimate = c(OR_coef2),
  StdError = c(OR_se2),
  CI_Lower = c(OR_CI_lb2),
  CI_Upper = c(OR_CI_ub2)
)

print(OR_results2)


##########################
# Questions 6 (b) 
##########################
expit <- function(x) {
  1 / (1 + exp(-x))
}

logit <- function(p) {
  log(p / (1 - p))
}

x <- seq(0.001, 0.999, length.out = 1000)
z <- expit(logit(x) + logistic_beta_lm)

plot(x, z, type = "l", col = "blue", lwd = 2,
     xlab = "E[Y|A=0,L]", ylab = "E[Y|A=1,L]",
     main = "expit(logit(E[Y|A=0,L]) + beta)")
abline(0, 1, col = "gray", lty = 2)  # identity line



