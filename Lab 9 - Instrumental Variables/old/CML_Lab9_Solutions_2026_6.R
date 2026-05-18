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
    log_lucent = log(lucent),
    log_igfbp3 = log(igfbp3)
  )

glimpse(data)

# ── Analysis variables ────────────────────────────────────────────────────────
Y <- data$log_lucent    # outcome
A <- data$log_igfbp3   # exposure
R <- data$igfbp3g_num  # instrument (numeric: AA=1, AC=2, CC=3)

# Numeric coding of genotype
data <- data %>%
  mutate(igfbp3g_num = as.numeric(igfbp3g))

R <- data$igfbp3g_num


################################################################################
##  Question 1 — Association of Lucent Area with Genotype
################################################################################

# ── Linear trend ──────────────────────────────────────────────────────────────
lm_q1 <- lm(log_lucent ~ igfbp3g_num, data = data)
summary(lm_q1)

# ── Factor model (no ordering assumed) ───────────────────────────────────────
lm_q1_cat <- lm(log_lucent ~ factor(igfbp3g), data = data)
summary(lm_q1_cat)

# ── Pearson correlation ───────────────────────────────────────────────────────
cor.test(data$log_lucent, data$igfbp3g_num)

# ── Visualise ────────────────────────────────────────────────────────────────
ggplot(data, aes(x = igfbp3g, y = log_lucent, fill = igfbp3g)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 16, width = 0.45) +
  scale_fill_manual(values = c("tomato", "steelblue", "forestgreen")) +
  theme_minimal() +
  labs(x = "Genotype (igfbp3g)", y = "log(lucent area)",
       title = "Log(Lucent Area) by IGFBP-3 Genotype") +
  theme(legend.position = "none")

# > Interpretation: The genotype is significantly associated with
# > log(lucent area) (p = 0.018). Each additional C allele corresponds to
# > a decrease of ~0.154 in log-lucent area. Under the IV assumptions,
# > this association propagates entirely through IGFBP-3 levels.


################################################################################
##  Question 2 — Distribution of IGFBP-3 and Relevance of the Instrument
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
lm_q2 <- lm(log_igfbp3 ~ igfbp3g_num, data = data)
summary(lm_q2)

lm_q2_cat <- lm(log_igfbp3 ~ factor(igfbp3g), data = data)
anova(lm_q2_cat)

# > Interpretation: The genotype is significantly associated with log(IGFBP-3)
# > (p = 0.009), supporting the relevance assumption. Each additional C allele
# > is associated with a decrease of ~0.063 in log(IGFBP-3). This is the only
# > IV assumption directly testable from data.


################################################################################
##  Question 3 — DAG and IV Assumptions
################################################################################

grViz("
digraph IV_DAG {
  graph [layout = neato, bgcolor = transparent, splines = curved]
  node [shape = plaintext, fontname = 'Helvetica', fontsize = 16]
  R [pos = '0,0!',   label = 'R\n(igfbp3g)']
  A [pos = '3,0!',   label = 'A\n(log igfbp3)']
  Y [pos = '6,0!',   label = 'Y\n(log lucent)']
  U [pos = '4.5,2!', label = 'U']
  L [pos = '3,-2!',  label = 'L\n(age, afbcat,\nbrfcat, parcat, ev_pil)']

  edge [arrowhead = vee, arrowsize = 0.8, penwidth = 1.6, color = '#444444']
  R -> A
  A -> Y
  U -> A
  U -> Y
  L -> A
  L -> Y

  edge [arrowhead = vee, arrowsize = 0.7, penwidth = 1.2,
        style = dashed, color = '#888888']
  R -> Y [style = dashed, constraint = false]
  R -> L [style = dashed]
}
")

# The three IV assumptions:
#
# 1. Relevance (R not-indep A | L):
#    The IGFBP3g genotype must be associated with serum IGFBP-3 after
#    adjusting for L. Verified empirically in Question 2 (p = 0.009).
#
# 2. Exclusion restriction (R indep Y | A, L, U):
#    The genotype should not affect breast tissue density except through
#    IGFBP-3. Potentially violated if the genotype influences parity or age
#    at first birth, which in turn affect breast density.
#
# 3. Randomization (R indep U | L):
#    Mendelian inheritance randomizes allele transmission, so the genotype
#    should be independent of unmeasured confounders after conditioning on L.
#    Not testable from data; relies on biological reasoning.


################################################################################
##  Question 4 — IV Estimate via EIF (No Covariates)
##
##  beta_q4_hat = sum_i (R_i - R_bar)(Y_i - Y_bar) /
##             sum_i (R_i - R_bar)(A_i - A_bar)
################################################################################

# ── Point estimate ────────────────────────────────────────────────────────────
beta_q4 <- mean((R - mean(R)) * (Y - mean(Y))) /
       mean((R - mean(R)) * (A - mean(A)))

# ── Standard error via influence function ─────────────────────────────────────
EIC_beta_q4 <- (R - mean(R)) * (Y - mean(Y) - beta_q4 * (A - mean(A))) /
           mean((R - mean(R)) * (A - mean(A)))
se_beta_q4 <- sd(EIC_beta_q4) / sqrt(N)

# ── 95% CI ────────────────────────────────────────────────────────────────────
beta_q4_lb <- beta_q4 - 1.96 * se_beta_q4
beta_q4_ub <- beta_q4 + 1.96 * se_beta_q4

# ── Results on geometric mean ratio scale (10% increase in A) ─────────────────
results_q4 <- tibble(
  Method   = "EIF (no covariates)",
  Estimate = round(exp(beta_q4    * log(1.1)), 3),
  CI_Lower = round(exp(beta_q4_lb * log(1.1)), 3),
  CI_Upper = round(exp(beta_q4_ub * log(1.1)), 3)
)

knitr::kable(results_q4, digits = 3,
             caption = "Q4 — IV estimate (no covariates): GMR for 10% increase in IGFBP-3")

# ── 2SLS comparison (no covariates) ──────────────────────────────────────────
iv_q4       <- ivreg(log_lucent ~ log_igfbp3 | igfbp3g_num, data = data)
iv_q4_ct    <- coeftest(iv_q4, vcov = sandwich)
iv_q4_ci    <- coefci(iv_q4, vcov = sandwich)["log_igfbp3", ]
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

# > Interpretation: A 10% increase in IGFBP-3 increases the geometric mean of
# > dense breast tissue area by ~26.1% [95% CI: 0.1% to 58.9%]. The CI is wide,
# > reflecting the small sample and modest instrument strength.


################################################################################
##  Question 5 — DML IV Estimate with Covariates (SuperLearner)
##
##  Confounder set: afbcat + parcat (blocks the R -> L pathway)
##  No cross-fitting given n = 139.
################################################################################

L_5 <- data %>% dplyr::select(afbcat, parcat)

SL.library <- c("SL.glm", "SL.ranger", "SL.earth")

# ── Nuisance models ───────────────────────────────────────────────────────────
sl_Y <- SuperLearner(Y = Y, X = L_5, SL.library = SL.library)
sl_A <- SuperLearner(Y = A, X = L_5, SL.library = SL.library)
sl_R <- SuperLearner(Y = R, X = L_5, SL.library = SL.library)

# ── Nuisance predictions ──────────────────────────────────────────────────────
m_Y <- predict(sl_Y, newdata = L_5, onlySL = TRUE)$pred
m_A <- predict(sl_A, newdata = L_5, onlySL = TRUE)$pred
m_R <- predict(sl_R, newdata = L_5, onlySL = TRUE)$pred

# ── EIF estimator ─────────────────────────────────────────────────────────────
beta_q5 <- mean((R - m_R) * (Y - m_Y)) /
         mean((R - m_R) * (A - m_A))

# ── Standard error ────────────────────────────────────────────────────────────
EIC_beta_q5 <- (R - m_R) * (Y - m_Y - beta_q5 * (A - m_A)) /
             mean((R - m_R) * (A - m_A))
se_beta_q5  <- sd(EIC_beta_q5) / sqrt(N)

# ── 95% CI ────────────────────────────────────────────────────────────────────
beta_q5_lb <- beta_q5 - 1.96 * se_beta_q5
beta_q5_ub <- beta_q5 + 1.96 * se_beta_q5

# ── Results ───────────────────────────────────────────────────────────────────
results_q5 <- tibble(
  Method   = "EIF (afbcat + parcat)",
  Estimate = round(exp(beta_q5    * log(1.1)), 3),
  CI_Lower = round(exp(beta_q5_lb * log(1.1)), 3),
  CI_Upper = round(exp(beta_q5_ub * log(1.1)), 3)
)

knitr::kable(results_q5, digits = 3,
             caption = "Q5 — EIF IV estimate with covariates: GMR for 10% increase in IGFBP-3")

# ── 2SLS comparison (afbcat + parcat) ────────────────────────────────────────
iv_q5         <- ivreg(log_lucent ~ log_igfbp3 + afbcat + parcat |
                         igfbp3g_num + afbcat + parcat, data = data)
iv_q5_ct      <- coeftest(iv_q5, vcov = sandwich)
iv_q5_ci      <- coefci(iv_q5, vcov = sandwich)["log_igfbp3", ]
beta_q5_2sls    <- iv_q5_ct["log_igfbp3", 1]
beta_q5_2sls_lb <- iv_q5_ci[1]
beta_q5_2sls_ub <- iv_q5_ci[2]

results_q5 <- rbind(results_q5, tibble(
  Method   = "2SLS (afbcat + parcat)",
  Estimate = round(exp(beta_q5_2sls    * log(1.1)), 3),
  CI_Lower = round(exp(beta_q5_2sls_lb * log(1.1)), 3),
  CI_Upper = round(exp(beta_q5_2sls_ub * log(1.1)), 3)
))

knitr::kable(results_q5, digits = 3,
             caption = "Q5 — EIF vs 2SLS (afbcat + parcat): GMR for 10% increase in IGFBP-3")

# > Interpretation: After conditioning on afbcat and parcat to block the
# > genotype -> confounders pathway, a 10% increase in IGFBP-3 increases
# > the geometric mean of lucent area by ~24.9% [95% CI: 3.4% to 50.9%].
# > The interval no longer includes zero.


################################################################################
##  Question 6 — Assumption-Lean Regression (Comparison)
##
##  Model: E[Y | A, L] = beta * A + g(L)  (partially linear)
##  Full confounder set L = (age, afbcat, brfcat, parcat, ev_pil)
##
##  Two approaches:
##   (a) EIF estimator
##   (b) OLS residual-on-residual with sandwich SE
################################################################################

L_6 <- data %>% dplyr::select(age, afbcat, brfcat, parcat, ev_pil)

# ── Nuisance models ───────────────────────────────────────────────────────────
sl_Y2 <- SuperLearner(Y = Y, X = L_6, SL.library = SL.library)
sl_A2 <- SuperLearner(Y = A, X = L_6, SL.library = SL.library)
sl_R2 <- SuperLearner(Y = R, X = L_6, SL.library = SL.library)

m_Y2 <- predict(sl_Y2, newdata = L_6, onlySL = TRUE)$pred
m_A2 <- predict(sl_A2, newdata = L_6, onlySL = TRUE)$pred
m_R2 <- predict(sl_R2, newdata = L_6, onlySL = TRUE)$pred

# ── EIF estimator (IV, full L) ────────────────────────────────────────────────
beta_q6 <- mean((R - m_R2) * (Y - m_Y2)) /
            mean((R - m_R2) * (A - m_A2))

EIC_beta_q6 <- (R - m_R2) * (Y - m_Y2 - beta_q6 * (A - m_A2)) /
            mean((R - m_R2) * (A - m_A2))
se_beta_q6   <- sd(EIC_beta_q6) / sqrt(N)

beta_q6_lb <- beta_q6 - 1.96 * se_beta_q6
beta_q6_ub <- beta_q6 + 1.96 * se_beta_q6


# ── Results ───────────────────────────────────────────────────────────────────
results_q6 <- tibble(
  Method   = "EIF (full L)",
  Estimate = round(exp(beta_q6    * log(1.1)), 3),
  CI_Lower = round(exp(beta_q6_lb * log(1.1)), 3),
  CI_Upper = round(exp(beta_q6_ub * log(1.1)), 3)
)

knitr::kable(results_q6, digits = 3,
             caption = "Q6 — Assumption-lean regression: GMR for 10% increase in IGFBP-3")

# ── 2SLS comparison (full L) ──────────────────────────────────────────────────
iv_q6        <- ivreg(log_lucent ~ log_igfbp3 + age + afbcat + brfcat + parcat + ev_pil |
                        igfbp3g_num + age + afbcat + brfcat + parcat + ev_pil, data = data)
iv_q6_ct     <- coeftest(iv_q6, vcov = sandwich)
iv_q6_ci     <- coefci(iv_q6, vcov = sandwich)["log_igfbp3", ]
beta_2sls    <- iv_q6_ct["log_igfbp3", 1]
beta_2sls_lb <- iv_q6_ci[1]
beta_2sls_ub <- iv_q6_ci[2]

results_q6 <- rbind(results_q6, tibble(
  Method   = "2SLS (full L)",
  Estimate = round(exp(beta_2sls    * log(1.1)), 3),
  CI_Lower = round(exp(beta_2sls_lb * log(1.1)), 3),
  CI_Upper = round(exp(beta_2sls_ub * log(1.1)), 3)
))

knitr::kable(results_q6, digits = 3,
             caption = "Q6 — EIF vs 2SLS (full L): GMR for 10% increase in IGFBP-3")

# > Interpretation: The assumption-lean regression estimates a much smaller
# > association — ~4.5% increase [95% CI: -0.2% to 9.3%], not significant.
# > Compared to the IV estimates (~25%), this suggests negative unmeasured
# > confounding that masks the causal effect in the observational analysis.


#############################################
##  Summary & Forest Plot
#############################################

results_all <- tibble(
  Method   = c(
    "No covariates — EIF",
    "No covariates — 2SLS",
    "afbcat + parcat — EIF",
    "afbcat + parcat — 2SLS",
    "Full L — EIF",
    "Full L — 2SLS"
  ),
  Type     = c("EIF", "2SLS", "EIF", "2SLS", "EIF", "2SLS"),
  Estimate = exp(c(beta_q4, beta_q4_2sls, beta_q5, beta_q5_2sls, beta_q6, beta_2sls) * log(1.1)),
  Lower    = exp(c(beta_q4_lb, beta_q4_2sls_lb, beta_q5_lb, beta_q5_2sls_lb, beta_q6_lb, beta_2sls_lb) * log(1.1)),
  Upper    = exp(c(beta_q4_ub, beta_q4_2sls_ub, beta_q5_ub, beta_q5_2sls_ub, beta_q6_ub, beta_2sls_ub) * log(1.1))
)

knitr::kable(
  results_all %>% mutate(across(where(is.numeric), ~round(.x, 3))),
  caption = "Summary — IV results"
)

# ── Forest plot ───────────────────────────────────────────────────────────────
results_all %>%
  mutate(Method = factor(Method, levels = rev(Method))) %>%
  ggplot(aes(x = Estimate, y = Method, color = Type)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0.25, linewidth = 0.8) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1, scale = 100)) +
  scale_color_manual(values = c("EIF" = "tomato", "2SLS" = "forestgreen", "Regression" = "steelblue")) +
  theme_minimal() +
  labs(
    x     = "Geometric mean ratio (10% increase in IGFBP-3)",
    y     = NULL,
    color = NULL,
    title = "IV Estimates — IGFBP-3 Effect on Breast Density"
  ) +
  theme(legend.position = "bottom")
