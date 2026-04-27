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

q.model <- 
p.model <- 

q <- 
p <- 

# --- LM with sandwich SE ---
model <- 
lm_result <- 

beta_lm <- 
se_lm <- 

beta_lm_lb <- 
beta_lm_ub <- 
  
# --- EIF ---
beta_eif <-
se_eif <- 

beta_eif_lb <-
beta_eif_ub <- 

# --- AIPW ---
library(npcausal) 
ate.aipw <- 

ate_aipw <- 
se_aipw <- 

ate_aipw_lb <- 
ate_aipw_ub <- 

# --- TMLE ---
library(tmle)
ate.tmle <- 

ate_tmle <- 
se_tmle <- 

ate_tmle_lb <- 
ate_tmle_ub <- 



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
# Try first with age
age <- L$age

age_summary <- 
se_age_summary <- 


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
    weights <- 
    
    # Weighted mean, SE and CIs
    weighted_mean <- 
    
    eif <- 
    weighted_se <- 
    
    weighted_CI_lb <- 
    weighted_CI_ub <- 
    
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

q.model2 <- 
p.model2 <- 

q2 <- 
p2 <- 


# --- LM with sandwich SE ---
model2 <- 
lm_result2 <- 

beta_lm2 <- 
se_lm2 <- 

beta_lm_lb2 <- 
beta_lm_ub2 <- 

# --- EIF ---
beta_eif2 <- 
se_eif2 <- 

beta_eif_lb2 <- 
beta_eif_ub2 <- 


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

q3I.model <- 
q3IC.model <- 

p3I.model <- 
p3IC.model <- 

q3I <- 
q3IC <- 

p3I <- 
p3IC <- 

q3 <- 
p3 <- 

# --- LM with sandwich SE ---
model3 <- 
lm_result3 <- 

beta_lm3 <- 
se_lm3 <- 

beta_lm_lb3 <- 
beta_lm_ub3 <- 

# --- EIF ---
beta_eif3 <- 
se_eif3 <- 

beta_eif_lb3 <- 
beta_eif_ub3 <- 

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

q4I.model <- 
q4IC.model <- 

p4I.model <- 
p4IC.model <- 

q4I <- 
q4IC <- 

p4I <- 
p4IC <- 

q4 <- 
p4 <- 

# --- LM with sandwich SE ---
model4 <- 
lm_result4 <- 

beta_lm4 <- 
se_lm4 <- 

beta_lm_lb4 <- 
beta_lm_ub4 <- 

# --- EIF ---
beta_eif4 <- 
se_eif4 <- 

beta_eif_lb4 <- 
beta_eif_ub4 <- 

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

Y_hat.model <- 
Y_hat <- 

#################
# Risk Ratio
#################

log_Y <- 

q_log.model <- 
q_log <- 

# --- LM with sandwich SE ---
pseudo_Y_log <- 

log_model <- 
log_lm_result <- 

log_beta_lm <- 
log_se_lm <- 

log_beta_lm_lb <- 
log_beta_lm_ub <- 

RR_coef <- 
RR_se <- 

RR_CI_lb <- 
RR_CI_ub <- 

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

logit_Y <- 

q_logit.model <- 
q_logit <- 

# --- LM with sandwich SE ---
pseudo_Y_logistic <- 

logistic_model <- 
logistic_lm_result <- 

logistic_beta_lm <- 
logistic_se_lm <- 

logistic_beta_lm_lb <- 
logistic_beta_lm_ub <- 
  
OR_coef <- 
OR_se <- 

OR_CI_lb <- 
OR_CI_ub <- 
  
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

y_hat.model <- 
y_hat <- 

#################
# Risk Ratio
#################

log_y <-

q2_log.model <-
q2_log <- 

# --- LM with sandwich SE ---
pseudo_y_log <- 

log_model2 <- 
log_lm_result2 <- 

log_beta_lm2 <- 
log_se_lm2 <- 

log_beta_lm_lb2 <- 
log_beta_lm_ub2 <- 

RR_coef2 <- 
RR_se2 <- 

RR_CI_lb2 <- 
RR_CI_ub2 <- 

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

logit_y <- 

q2_logit.model <- 
q2_logit <- 

# --- LM with sandwich SE ---
pseudo_y_logistic <- 

logistic_model2 <- 
logistic_lm_result2 <- 

logistic_beta_lm2 <- 
logistic_se_lm2 <- 

logistic_beta_lm_lb2 <- 
logistic_beta_lm_ub2 <- 

OR_coef2 <- 
OR_se2 <- 

OR_CI_lb2 <- 
OR_CI_ub2 <- 

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



