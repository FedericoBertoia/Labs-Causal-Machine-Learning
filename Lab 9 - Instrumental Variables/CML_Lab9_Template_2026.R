#############################################
##  Instrumental Variable Analysis
##  IGFBP-3 and Dense Breast Tissue (n = 139)
#############################################

# ── Packages ───────────────────────────────────────────────────────────────────
library(foreign)
library(dplyr)
library(ggplot2)
library(patchwork)
library(tibble)
library(SuperLearner)
library(lmtest)
library(sandwich)
library(DiagrammeR)
library(AER)

set.seed(123)


#############################################
##  Data Loading & Preparation
#############################################

data <- read.dta("igf.dta")
N    <- nrow(data)

# ── Log-transform exposure and outcome ────────────────────────────────────────
data <- data %>%
  mutate(
    log_lucent  = log(lucent),
    log_igfbp3  = log(igfbp3),
    igfbp3g_num = as.numeric(igfbp3g)   # AA=1, AC=2, CC=3
  )

glimpse(data)

# ── Analysis variables ────────────────────────────────────────────────────────
Y <- data$log_lucent    # outcome
A <- data$log_igfbp3   # exposure
R <- data$igfbp3g_num  # instrument (numeric)


################################################################################
##  Question 1 — Association of Lucent Area with Genotype
##
##  Is the instrument (genotype) associated with the outcome (log lucent area)?
##  Fit a linear model and visualise.
################################################################################

# ── Linear trend ──────────────────────────────────────────────────────────────
lm_q1 <- lm(log_lucent ~ ..., data = data)
summary(lm_q1)

# ── Factor model (no ordering assumed) ───────────────────────────────────────
lm_q1_cat <- lm(log_lucent ~ ..., data = data)
summary(lm_q1_cat)

# ── Pearson correlation ───────────────────────────────────────────────────────
cor.test(...)

# ── Visualise ────────────────────────────────────────────────────────────────
ggplot(data, aes(x = igfbp3g, y = log_lucent, fill = igfbp3g)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 16, width = 0.45) +
  scale_fill_manual(values = c("tomato", "steelblue", "forestgreen")) +
  theme_minimal() +
  labs(x = "Genotype (igfbp3g)", y = "log(lucent area)",
       title = "Log(Lucent Area) by IGFBP-3 Genotype") +
  theme(legend.position = "none")

# Interpret the association. What role does this association play under
# the IV assumptions?


################################################################################
##  Question 2 — Distribution of IGFBP-3 and Relevance of the Instrument
##
##  Test whether the genotype (R) is associated with the exposure (log IGFBP-3).
##  This is the only IV assumption directly testable from data.
################################################################################

# ── 2.1: Histograms (raw and log-scale) ──────────────────────────────────────
p1 <- ggplot(data, aes(x = igfbp3)) +
  geom_histogram(bins = 10, fill = "steelblue", alpha = 0.7, color = "white") +
  theme_minimal() +
  labs(x = "igfbp3", y = "Frequency", title = "Histogram of igfbp3")

p2 <- ggplot(data, aes(x = log_igfbp3)) +
  geom_histogram(bins = 10, fill = "steelblue", alpha = 0.7, color = "white") +
  theme_minimal() +
  labs(x = "log_igfbp3", y = "Frequency", title = "Histogram of log(igfbp3)")

p1 + p2

# ── 2.2: Distribution by genotype ─────────────────────────────────────────────
data %>%
  group_by(igfbp3g) %>%
  summarise(
    n    = n(),
    mean = round(mean(log_igfbp3), 3),
    sd   = round(sd(log_igfbp3),   3),
    min  = round(min(log_igfbp3),  3),
    max  = round(max(log_igfbp3),  3)
  )

ggplot(data, aes(x = log_igfbp3, fill = igfbp3g, colour = igfbp3g)) +
  geom_histogram(alpha = 0.4, position = "identity", bins = 10) +
  scale_fill_manual(values   = c("tomato", "steelblue", "forestgreen")) +
  scale_colour_manual(values = c("tomato", "steelblue", "forestgreen")) +
  theme_minimal() +
  labs(x = "log_igfbp3", y = "Count",
       title = "Distribution of log(IGFBP-3) by Genotype",
       fill = "igfbp3g", colour = "igfbp3g")

# ── 2.3: Testing relevance ────────────────────────────────────────────────────
lm_q2 <- lm(log_igfbp3 ~ ..., data = data)
summary(lm_q2)

lm_q2_cat <- lm(log_igfbp3 ~ ..., data = data)
anova(lm_q2_cat)

# Interpret the results. Which IV assumption does this test address?
# Is this assumption sufficient to establish causality on its own?


################################################################################
##  Question 3 — DAG and IV Assumptions
##
##  Draw a DAG for this study and state the three IV assumptions in the
##  context of the IGFBP-3 / breast density analysis.
################################################################################

# ── DAG ───────────────────────────────────────────────────────────────────────
grViz("
digraph IV_DAG {
  graph [layout = neato, bgcolor = transparent, splines = curved]
  node [shape = plaintext, fontname = 'Helvetica', fontsize = 16]

  # Nodes — fill in positions and labels
  R [pos = '0,0!',   label = 'R\n(igfbp3g)']
  A [pos = '3,0!',   label = 'A\n(log igfbp3)']
  Y [pos = '6,0!',   label = 'Y\n(log lucent)']
  U [pos = '4.5,2!', label = 'U']
  L [pos = '3,-2!',  label = 'L\n(...)']

  edge [arrowhead = vee, arrowsize = 0.8, penwidth = 1.6, color = '#444444']
  # Add edges for causal paths
  ...
}
")

# State and discuss the three IV assumptions:
#
# 1. Relevance (R not-indep A | L):
#    ...
#
# 2. Exclusion restriction (R indep Y | A, L, U):
#    ...
#
# 3. Randomization / independence (R indep U | L):
#    ...
#
# Which of these can be empirically tested? Which rely on subject-matter knowledge?


################################################################################
##  Question 4 — IV Estimate via EIF (No Covariates)
################################################################################

# ── Point estimate ────────────────────────────────────────────────────────────
beta_q4 <- ... / ...

# ── Standard error via influence function ────────────────────────────────────
EIC_beta_q4 <- ...
se_beta_q4  <- sd(EIC_beta_q4) / sqrt(N)

# ── 95% CI ────────────────────────────────────────────────────────────────────
beta_q4_lb <- beta_q4 - 1.96 * se_beta_q4
beta_q4_ub <- beta_q4 + 1.96 * se_beta_q4

# ── Results on geometric mean ratio scale (10% increase in A) ─────────────────
# Hint: GMR for a 10% increase = exp(beta * log(1.1))
results_q4 <- tibble(
  Method   = "EIF (no covariates)",
  Estimate = round(exp(beta_q4    * log(1.1)), 3),
  CI_Lower = round(exp(beta_q4_lb * log(1.1)), 3),
  CI_Upper = round(exp(beta_q4_ub * log(1.1)), 3)
)

# ── 2SLS comparison (no covariates) ──────────────────────────────────────────
iv_q4    <- ivreg(log_lucent ~ log_igfbp3 | igfbp3g_num, data = data)
iv_q4_ct <- coeftest(iv_q4, vcov = sandwich)
iv_q4_ci <- coefci(iv_q4,   vcov = sandwich)["log_igfbp3", ]

beta_q4_2sls    <- iv_q4_ct["log_igfbp3", 1]
beta_q4_2sls_lb <- iv_q4_ci[1]
beta_q4_2sls_ub <- iv_q4_ci[2]

results_q4 <- rbind(results_q4, tibble(
  Method   = "2SLS (no covariates)",
  Estimate = round(exp(beta_q4_2sls    * log(1.1)), 3),
  CI_Lower = round(exp(beta_q4_2sls_lb * log(1.1)), 3),
  CI_Upper = round(exp(beta_q4_2sls_ub * log(1.1)), 3)
))

knitr::kable(results_q4, digits = 3,
             caption = "Q4 — EIF vs 2SLS (no covariates): GMR for 10% increase in IGFBP-3")

# Interpret the estimate and the width of the CI.
# Does the EIF estimator agree with 2SLS? Why or why not?


################################################################################
##  Question 5 — DML IV Estimate with Covariates (SuperLearner)
################################################################################

L_5 <- data %>% dplyr::select(...)

SL.library <- c("SL.glm", "SL.ranger", "SL.earth")

# ── Nuisance models ───────────────────────────────────────────────────────────
sl_Y_q5 <- SuperLearner(Y = ..., X = L_5, SL.library = SL.library)
sl_A_q5 <- SuperLearner(Y = ..., X = L_5, SL.library = SL.library)
sl_R_q5 <- SuperLearner(Y = ..., X = L_5, SL.library = SL.library)

# ── Nuisance predictions ──────────────────────────────────────────────────────
m_Y_q5 <- predict(sl_Y_q5, newdata = L_5, onlySL = TRUE)$pred
m_A_q5 <- predict(sl_A_q5, newdata = L_5, onlySL = TRUE)$pred
m_R_q5 <- predict(sl_R_q5, newdata = L_5, onlySL = TRUE)$pred

# ── EIF estimator ─────────────────────────────────────────────────────────────
beta_q5 <- ... / ...

# ── Standard error ────────────────────────────────────────────────────────────
EIC_beta_q5 <- ...
se_beta_q5  <- sd(EIC_beta_q5) / sqrt(N)

# ── 95% CI ────────────────────────────────────────────────────────────────────
beta_q5_lb <- beta_q5 - 1.96 * se_beta_q5
beta_q5_ub <- beta_q5 + 1.96 * se_beta_q5

# ── Results ───────────────────────────────────────────────────────────────────
results_q5 <- tibble(
  Method   = "EIF (covariates)",
  Estimate = round(exp(beta_q5    * log(1.1)), 3),
  CI_Lower = round(exp(beta_q5_lb * log(1.1)), 3),
  CI_Upper = round(exp(beta_q5_ub * log(1.1)), 3)
)

# ── 2SLS comparison (same covariates) ────────────────────────────────────────
iv_q5    <- ivreg(log_lucent ~ log_igfbp3 + ... | igfbp3g_num + ..., data = data)
iv_q5_ct <- coeftest(iv_q5, vcov = sandwich)
iv_q5_ci <- coefci(iv_q5,   vcov = sandwich)["log_igfbp3", ]

beta_q5_2sls    <- iv_q5_ct["log_igfbp3", 1]
beta_q5_2sls_lb <- iv_q5_ci[1]
beta_q5_2sls_ub <- iv_q5_ci[2]

results_q5 <- rbind(results_q5, tibble(
  Method   = "2SLS (covariates)",
  Estimate = round(exp(beta_q5_2sls    * log(1.1)), 3),
  CI_Lower = round(exp(beta_q5_2sls_lb * log(1.1)), 3),
  CI_Upper = round(exp(beta_q5_2sls_ub * log(1.1)), 3)
))

knitr::kable(results_q5, digits = 3,
             caption = "Q5 — EIF vs 2SLS (with covariates): GMR for 10% increase in IGFBP-3")

# How does adjusting for covariates change the estimate and the CI?


################################################################################
##  Question 6 — Assumption-Lean Regression (No Instrument)
################################################################################

L_6 <- data %>% dplyr::select(age, afbcat, brfcat, parcat, ev_pil)

# ── Nuisance models (Y and A only — no instrument) ───────────────────────────
sl_Y_q6 <- SuperLearner(Y = ..., X = L_6, SL.library = SL.library)
sl_A_q6 <- SuperLearner(Y = ..., X = L_6, SL.library = SL.library)

m_Y_q6 <- predict(sl_Y_q6, newdata = L_6, onlySL = TRUE)$pred
m_A_q6 <- predict(sl_A_q6, newdata = L_6, onlySL = TRUE)$pred

# ── EIF estimator ──────────────────────────────────────────────────
A_res <- A - m_A_q6
Y_res <- Y - m_Y_q6

beta_q6 <- ... / ...

# ── Standard error via influence function ────────────────────────────────────
EIC_beta_q6 <- ...
se_beta_q6  <- sd(EIC_beta_q6) / sqrt(N)

# ── 95% CI ────────────────────────────────────────────────────────────────────
beta_q6_lb <- beta_q6 - 1.96 * se_beta_q6
beta_q6_ub <- beta_q6 + 1.96 * se_beta_q6

# ── Results ───────────────────────────────────────────────────────────────────
results_q6 <- tibble(
  Method   = "Assumption-lean regression (full L)",
  Estimate = round(exp(beta_q6    * log(1.1)), 3),
  CI_Lower = round(exp(beta_q6_lb * log(1.1)), 3),
  CI_Upper = round(exp(beta_q6_ub * log(1.1)), 3)
)

# Compare the regression estimate to the IV estimates from Q4 and Q5.
# What does the discrepancy suggest about unmeasured confounding?
# Under what conditions would the regression and IV estimates agree?


################################################################################
##  Forest Plot — All Estimates
################################################################################

results_all <- tibble(
  Method = c(
    "IV — EIF (no covariates)",
    "IV — 2SLS (no covariates)",
    "IV — EIF (covariates)",
    "IV — 2SLS (covariates)",
    "Regression — OLS resid-on-resid (full L)"
  ),
  Type     = c("IV", "IV", "IV", "IV", "Regression"),
  Estimate = exp(c(beta_q4, beta_q4_2sls, beta_q5, beta_q5_2sls, beta_q6) * log(1.1)),
  Lower    = exp(c(beta_q4_lb, beta_q4_2sls_lb, beta_q5_lb, beta_q5_2sls_lb, beta_q6_lb) * log(1.1)),
  Upper    = exp(c(beta_q4_ub, beta_q4_2sls_ub, beta_q5_ub, beta_q5_2sls_ub, beta_q6_ub) * log(1.1))
)

results_all %>%
  mutate(Method = factor(Method, levels = rev(Method))) %>%
  ggplot(aes(x = Estimate, y = Method, colour = Type)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey60") +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0.25, linewidth = 0.8) +
  geom_point(size = 3) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1, scale = 100)) +
  scale_colour_manual(values = c("IV" = "tomato", "Regression" = "steelblue")) +
  theme_minimal(base_size = 13) +
  labs(
    x      = "Geometric mean ratio (10% increase in IGFBP-3)",
    y      = NULL,
    colour = NULL,
    title  = "IV vs Regression Estimates — IGFBP-3 Effect on Breast Density"
  ) +
  theme(legend.position = "bottom")
