library(foreign)

# Load the data

file_path <-
  
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
anova(lm_model_cat)

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
Y <- 

# The treatment variable
A <- 

# The instrumental variable
R <- 

# Compute the estimator psi 
psi <- 

# Compute standard errors for psi using the efficient influence curves
EIC.psi <- 
se.psi <- 

# Compute the confidence intervals for psi
psi.CI.low <- 
psi.CI.high <- 



results <- tibble::tibble(
  Method = c("Eff. Influence Function"),
  Estimate = format(round(exp(psi * log(1.1)), 3), nsmall = 3),
  CI_Lower = format(round(exp(psi.CI.low * log(1.1)), 3), nsmall = 3),
  CI_Upper = format(round(exp(psi.CI.high * log(1.1)), 3), nsmall = 3)
)

print(results)


##################################
### Question 5
##################################

library(SuperLearner)

SL.library <-  

# Set of covariates assumed to mediate the effect of the instrument on the outcome
L <- 

# Fit nuisance parameters, i.e. prediction models (given covariates L) for the outcome E[Y|L], for the exposure E[A|L] and for the instrument E[R|L]
YL.sl <- 
AL.sl <- 
RL.sl <- 

# Make predictions using fitted models
L.YL = 
L.AL = 
L.RL = 

# Compute the estimator psi with covariates
psi.L <- 

# Compute standard errors for psi.L using the efficient influence curve
EIC.psi.L <- 
se.psi.L <- 

# Compute the confidence intervals for psi.L
psi.L.CI.low <- 
psi.L.CI.high <- 



results2 <- tibble::tibble(
  Method = c("Eff. Influence Function"),
  Estimate = format(round(exp(psi.L * log(1.1)), 3), nsmall = 3),
  CI_Lower = format(round(exp(psi.L.CI.low * log(1.1)), 3), nsmall = 3),
  CI_Upper = format(round(exp(psi.L.CI.high * log(1.1)), 3), nsmall = 3)
)

print(results2)

##################################
### Question 6
##################################

library(lmtest)
library(sandwich)

# Assumption-Lean Linear Regression (
# We are interested in estimating the association between the outcome (i.e. log(lucent)) 
# and the exposure of interest (i.e. igfbp3).

# Set of confounders 
L2 <- 

# Fit nuisance parameters, i.e. the outcome prediction model E[Y|L] and the propensity score model E[A|L] 
YL2.sl <- 
AL2.sl <- 

# Make predictions using fitted models
YL2 = 
AL2 = 

# --- EIF ---

beta <- 

EIC.beta <- 
se.beta <- 

beta.CI.low <- 
beta.CI.high <- 

# --- LM with sandwich SE ---
lm_model <- 
lm_result <- 

beta_lm<- 
se_lm <- 

beta_lm.CI.low <- 
beta_lm.CI.high <- 




results3 <- tibble::tibble(
  Method = c("Linear Model", "Eff. Influence Function"),
  Estimate = format(round(exp(c(beta_lm, beta) * log(1.1)), 3), nsmall = 3),
  CI_Lower = format(round(exp(c(beta_lm.CI.low, beta.CI.low) * log(1.1)), 3), nsmall = 3),
  CI_Upper = format(round(exp(c(beta_lm.CI.high, beta.CI.high) * log(1.1)), 3), nsmall = 3)
)

print(results3)




