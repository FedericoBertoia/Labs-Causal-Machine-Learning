library(foreign)

# Load the data

file_path <- "C:/Users/fbertoia/OneDrive - UGent/Desktop/PhD Statistical Data Analysis/1. Causal Machine Learning - Stijn Vanstelaandt/Labs/Lab 9 - Instrumental Variables/igf.dta"
file_path2 <- "C:/Users/feder/OneDrive - UGent/Desktop/PhD Statistical Data Analysis/1. Causal Machine Learning - Stijn Vanstelaandt/Labs/Lab 9 - Instrumental Variables/igf.dta"

data <- read.dta(file_path)
attach(data)
N <- dim(data)[1]

set.seed(123)

############################################################################

### Data description
# =========================================================================
#
# id          Patient's ID    
# igfbp3g     Putative functional gene for the protein IGFBP-3 (with genotypes "A/A", "A/C" and "C/C").
# age         Patient's age.
# brfcat      Breast feeding duration
# ev_pil      Indicator whether ever taken oral contraceptives. 
# ev_smok     Indicator whether ever smoked
# bmi         Patient's BMI
# igfbp3      Serum IGFBP-3 levels
# lucent      Total area of dense tissue in the breast
# agecat      Patient's age (categorized)
# parcat      Number of children
# afbcat      Age at fist birth (categorized); nulliparous - never having given birth before
# bmicat      Patient's BMI (categorized)

############################################################################

##################################
### Question 1
##################################

# Association between (log) lucent area and genotype
data$log_lucent <- log(lucent)
attach(data)

cor(log_lucent, as.numeric(igfbp3g))
cor.test(log_lucent, as.numeric(igfbp3g), method=c("pearson")) 

model <- lm(log_lucent ~ as.numeric(igfbp3g), data = data)
summary(model)

# As factor
lm_model_cat <- lm(log_lucent ~ factor(igfbp3g), data = data)
summary(lm_model_cat)

##################################
### Question 2
##################################

library(purrr)
library(ggplot2)

# Histogram of log(serum IGFBP-3)
data$log_igfbp3 <- log(igfbp3)
attach(data)

hist(igfbp3)
hist(log_igfbp3)

# Summary the distribution of log(serum IGFBP-3) by genotype
data %>% 
  split(.$igfbp3g) %>% 
  map(summary)


# Histogram of log(serum IGFBP-3) by genotype
ggplot(data, aes(x=log_igfbp3, fill=igfbp3g, colour = igfbp3g)) +
  geom_histogram(alpha=0.4, position = 'identity') 

# Test relevance of the instrument

model2 <- lm(log_igfbp3 ~ as.numeric(igfbp3g), data = data)
summary(model2)

# As factor
lm_model_cat2 <- lm(log_igfbp3 ~ factor(igfbp3g), data = data)
summary(lm_model_cat2)
anova(lm_model_cat2)

##################################
### Question 4
##################################

# The outcome variable 
Y <- log_lucent

# The treatment variable
A <- log_igfbp3

# The instrumental variable
R <- as.numeric(igfbp3g)

# Compute the estimator psi 
psi <- mean((R - mean(R))*(Y - mean(Y))) / mean((R - mean(R))*(A - mean(A)))

# Compute standard errors for psi using the efficient influence curves
EIC.psi <- (R - mean(R))*(Y - mean(Y) - psi*(A - mean(A))) / mean((R - mean(R))*(A - mean(A)))
se.psi <- (1/sqrt(N)) * sd(EIC.psi)

# Compute the confidence intervals for psi
psi.CI.low <- psi - 1.96 * se.psi
psi.CI.high <- psi + 1.96 * se.psi



results <- tibble::tibble(
  Method = c("Eff. Influence Function"),
  Estimate = format(round(exp(psi * log(1.1)), 3), nsmall = 3),
  CI_Lower = format(round(exp(psi.CI.low * log(1.1)), 3), nsmall = 3),
  CI_Upper = format(round(exp(psi.CI.high * log(1.1)), 3), nsmall = 3)
)

print(results)



# Interpretation: A 10% increase in level of serum IGFBP-3 increases the geometric mean 
# of the lucent area by 26.15% (= exp(2.44 * log(1.1))) [95% CI: 0.11% to 58.94%].

##################################
### Question 5
##################################

library(SuperLearner)

SL.library <-  c("SL.glm", "SL.randomForest")

# Set of covariates assumed to mediate the effect of the instrument on the outcome
L <- data[,c("afbcat", "parcat")]

# Fit nuisance parameters, i.e. prediction models (given covariates L) for the outcome E[Y|L], for the exposure E[A|L] and for the instrument E[R|L]
YL.sl <- SuperLearner(Y = Y, X = L, SL.library = SL.library)
AL.sl <- SuperLearner(Y = A, X = L, SL.library = SL.library)
RL.sl <- SuperLearner(Y = R, X = L, SL.library = SL.library)

# Make predictions using fitted models
L.YL = predict(YL.sl, newdata = L, onlySL = TRUE)$pred
L.AL = predict(AL.sl, newdata = L)$pred
L.RL = predict(RL.sl, newdata = L)$pred

# Compute the estimator psi with covariates
psi.L <- mean((R - L.RL)*(Y - L.YL)) / mean((R - L.RL)*(A - L.AL)); 

# Compute standard errors for psi.L using the efficient influence curve
EIC.psi.L <- (R - L.RL)*(Y - L.YL - psi.L*(A - L.AL)) / mean((R - L.RL)*(A - L.AL))
se.psi.L <- (1/sqrt(N)) * sd(EIC.psi.L); 

# Compute the confidence intervals for psi.L
psi.L.CI.low <- psi.L - 1.96 * se.psi.L
psi.L.CI.high <- psi.L + 1.96 * se.psi.L



results2 <- tibble::tibble(
  Method = c("Eff. Influence Function"),
  Estimate = format(round(exp(psi.L * log(1.1)), 3), nsmall = 3),
  CI_Lower = format(round(exp(psi.L.CI.low * log(1.1)), 3), nsmall = 3),
  CI_Upper = format(round(exp(psi.L.CI.high * log(1.1)), 3), nsmall = 3)
)

print(results2)

# Interpretation: A 10% increase in level of serum IGFBP-3 increases the geometric mean 
# of the lucent area by 24.9% (= exp(2.33 * log(1.1))) [95% CI: 3.4% to 50.9%].

##################################
### Question 6
##################################

library(lmtest)
library(sandwich)

# Assumption-Lean Linear Regression (
# We are interested in estimating the association between the outcome (i.e. log(lucent)) 
# and the exposure of interest (i.e. igfbp3).

# Set of confounders 
L2 <- data[,c("age", "afbcat", "brfcat", "parcat", "ev_pil")]

# Fit nuisance parameters, i.e. the outcome prediction model E[Y|L] and the propensity score model E[A|L] 
YL2.sl <- SuperLearner(Y = Y, X = L2, SL.library = SL.library)
AL2.sl <- SuperLearner(Y = A, X = L2, SL.library = SL.library)

# Make predictions using fitted models
YL2 = predict(YL2.sl, newdata = L2)$pred
AL2 = predict(AL2.sl, newdata = L2)$pred

# --- EIF ---

beta <- mean((A-AL2)*(Y-YL2)) / mean((A-AL2)^2)

EIC.beta <- (A - AL2)*(Y - YL2 - beta*(A - AL2)) / mean((A - AL2)^2)
se.beta <- (1/sqrt(N)) * sd(EIC.beta)

beta.CI.low <- beta - 1.96 * se.beta
beta.CI.high <- beta + 1.96 * se.beta

# --- LM with sandwich SE ---
lm_model <- lm(I(Y-YL2) ~ -1 + I(A-AL2))
lm_result <- coeftest(lm_model, vcov = sandwich)

beta_lm<- lm_result[1, 1]
se_lm <- lm_result[1, 2]

beta_lm.CI.low <- beta_lm - 1.96 * se_lm
beta_lm.CI.high <- beta_lm + 1.96 * se_lm




results3 <- tibble::tibble(
  Method = c("Linear Model", "Eff. Influence Function"),
  Estimate = format(round(exp(c(beta_lm, beta) * log(1.1)), 3), nsmall = 3),
  CI_Lower = format(round(exp(c(beta_lm.CI.low, beta.CI.low) * log(1.1)), 3), nsmall = 3),
  CI_Upper = format(round(exp(c(beta_lm.CI.high, beta.CI.high) * log(1.1)), 3), nsmall = 3)
)

print(results3)


# Interpretation: A 10% increase in level of serum IGFBP-3 increases the geometric mean 
# of the lucent area by 4.5% (= exp(0.46 * log(1.1))) [95% CI: -0.2% to 9.5%].


